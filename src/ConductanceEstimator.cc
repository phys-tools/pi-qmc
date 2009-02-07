// $Id$
/*  Copyright (C) 2004-2008 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ConductanceEstimator.h"
#include "SimulationInfo.h"
#include "SuperCell.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include <fstream>

ConductanceEstimator::ConductanceEstimator(const SimulationInfo& simInfo, 
    const int nfreq, Species* species, bool useSpeciesTensor,
    int idim, bool useCharge, MPIManager *mpi, int norder)
  : BlitzArrayBlkdEst<6>("conductance", 
      IVecN(norder, (idim==-1)?NDIM:1, (idim==-1)?NDIM:1,
            useSpeciesTensor?simInfo.getNSpecies():1,
            useSpeciesTensor?simInfo.getNSpecies():1,nfreq),true),
    npart(species?species->count:simInfo.getNPart()),
    ifirst(species?species->ifirst:0), nslice(simInfo.getNSlice()), 
    nfreq(nfreq), tau(simInfo.getTau()), idim(idim),
    temp(value.shape()+blitz::TinyVector<int,6>(0,0,0,0,0,nslice-nfreq)),
    buff(value.shape()+blitz::TinyVector<int,6>(0,0,0,0,0,nslice-nfreq)),
    mpi(mpi), charge(npart), ispecies(npart), cell(*simInfo.getSuperCell()),
    norder(norder) {
  // Set up charge array (or set to ones if not charge current).
  if (useCharge) {
    for (int i=0; i<npart; ++i) charge(i)=simInfo.getPartSpecies(i).charge;
  } else {
    charge=1.;
  }
  // Set up species index if making a species tensor.
  if (useSpeciesTensor) {
    const Species *species=&simInfo.getPartSpecies(ifirst);
    int label=0;
    for (int i=0; i<npart; ++i) {
      const Species *temp=&simInfo.getPartSpecies(i+ifirst);
      if (temp!=species) ++label;
      ispecies(i)=label;
      species=temp;
    }
  } else {
    ispecies=0;
  }
  // Set up FFT's.
  fftw_complex *ptr = (fftw_complex*)temp.data();
  int howmany=value.size()/nfreq/norder;// fft each component
  int istride=1; // time or freq is stored contiguous
  fwd = fftw_plan_many_dft(1,&nslice,howmany,
                           ptr,0,istride,nslice,
                           ptr,0,istride,nslice,
                           FFTW_FORWARD,FFTW_MEASURE);
  rev = fftw_plan_many_dft(1,&nslice,howmany,
                           ptr,0,istride,nslice,
                           ptr,0,istride,nslice,
                           FFTW_BACKWARD,FFTW_MEASURE);
}

ConductanceEstimator::~ConductanceEstimator() {
  fftw_destroy_plan(rev);
  fftw_destroy_plan(fwd);
}

void ConductanceEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0; n=0;
}

void ConductanceEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec diff=end-start;
    if (idim==-1) {
      for (int i=0; i<NDIM; ++i) {
        temp(0,i,i,ispecies(ipart-ifirst),ispecies(ipart-ifirst),islice)
          +=charge(ipart)*cell.pbc(diff)[i];
      }
    } else {
      temp(0,0,0,ispecies(ipart-ifirst),ispecies(ipart-ifirst),islice)
        +=charge(ipart)*cell.pbc(diff)[idim];
    }
  }
}

void ConductanceEstimator::endCalc(const int lnslice) {
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    buff=0;
    mpi->getWorkerComm().Reduce(temp.data(),buff.data(),
                                2*temp.size(),MPI::DOUBLE,MPI::SUM,0);
    if (workerID==0) temp=buff;
  }
#endif
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    blitz::Range freqs(0,nfreq-1), all(blitz::Range::all()), zero(0,0);
    fftw_execute(fwd);
    int ndim=value.extent(1), nspecies=value.extent(3);
    // First handle higher orders.
    for (int iorder=2; iorder<norder+1; ++iorder) { 
      int nnfreq = nslice/2/iorder;
      if (nnfreq>nfreq) nnfreq=nfreq;
      for (int id=0; id<ndim; ++id) { 
        for (int jd=0; jd<ndim; ++jd) { 
          for (int i=0; i<nspecies; ++i) { 
            for (int j=0; j<nspecies; ++j) { 
              if (iorder%2==0) {
                for (int ifreq=0; ifreq<nnfreq; ++ifreq) {
                  temp(iorder-1,id,jd,i,j,ifreq) = imag(
                     pow(temp(0,jd,jd,j,j,ifreq),iorder)
                   *conj(temp(0,id,id,i,i,iorder*ifreq)));
                }
              } else {
                for (int ifreq=0; ifreq<nnfreq; ++ifreq) {
                  temp(iorder-1,id,jd,i,j,ifreq) = real(
                     pow(temp(0,jd,jd,j,j,ifreq),iorder)
                   *conj(temp(0,id,id,i,i,iorder*ifreq)));
                }
              }
            }
          }
        }
      }
    }
    // Then handle order 1 (more difficult because data is in place.
    // First calculate elements between species.
    for (int ispec=0; ispec<nspecies; ++ispec) {
      blitz::Range ispc(ispec,ispec);
      for (int jspec=0; jspec<nspecies; ++jspec) {
        blitz::Range jspc(jspec,jspec);
        if (ispec!=jspec) {
          if (idim==-1) {
            for (int i=0; i<NDIM; ++i) {
              blitz::Range idim=blitz::Range(i,i);
              for (int j=0; j<NDIM; ++j) {
                blitz::Range jdim=blitz::Range(j,j);
                temp(zero,idim,jdim,ispc,jspc,all)
                       =temp(zero,idim,idim,ispc,ispc,all)
                  *conj(temp(zero,jdim,jdim,jspc,jspc,all));
              }
            }
          } else {
            temp(zero,zero,zero,ispc,jspc,all)
                   =temp(zero,zero,zero,ispc,ispc,all)
              *conj(temp(zero,zero,zero,jspc,jspc,all));
          }
        }
      }
    }
    // Then handle elements within a species.
    for (int ispec=0; ispec<nspecies; ++ispec) {
      blitz::Range spec(ispec,ispec);
      if (idim==-1) {
        for (int i=0; i<NDIM; ++i) {
          blitz::Range idim(i,i);
          for (int j=i+1; j<NDIM; ++j) {
            blitz::Range jdim(j,j);
            temp(zero,idim,jdim,spec,spec,all)
                  = temp(zero,idim,idim,spec,spec,all)
              *conj(temp(zero,jdim,jdim,spec,spec,all));
            temp(zero,jdim,idim,spec,spec,all)
              = conj(temp(zero,idim,jdim,spec,spec,all));
          }
        }
        for (int i=0; i<NDIM; ++i) {
          blitz::Range idim(i,i);
          temp(zero,idim,idim,spec,spec,all)
            *= conj(temp(zero,idim,idim,spec,spec,all));
        }
      } else {
          temp(zero,all,all,spec,spec,all)
            *=conj(temp(zero,all,all,spec,spec,all));
      }
    }
    double betainv=1./(tau*nslice);
    value -= real(temp(all,all,all,all,all,freqs))*betainv;
    norm+=1;
  }
}
