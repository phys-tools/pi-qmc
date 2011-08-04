// $Id$
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#include "VIndEstimator.h"
#include "SimulationInfo.h"
#include "CoulombAction.h"
#include "SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "stats/MPIManager.h"

VIndEstimator::VIndEstimator(const SimulationInfo& simInfo,
  const CoulombAction* coulAction, const int nfreq, 
  const int nbin, const int nstride, MPIManager *mpi)
  : BlitzArrayBlkdEst<2>("vind","dynamic-array/induced-voltage",
                         IVecN(nbin,nfreq),true), 
    coulAction(coulAction), npart(simInfo.getNPart()), 
    nslice(simInfo.getNSlice()), nfreq(nfreq), nbin(nbin), 
    nstride(nstride), tau(simInfo.getTau()), tauinv(1./tau), 
    massinv(1./simInfo.getSpecies(0).mass),
    dx(simInfo.getSuperCell()->a[0]/nbin), dxinv(1/dx),q(npart),
    temp(nbin,nslice/nstride), 
    buff(nbin,nslice/nstride), 
    temp2(nslice/nstride), mpi(mpi) {
  fftw_complex *ptr = (fftw_complex*)temp.data();
  fftw_complex *ptr2 = (fftw_complex*)temp2.data();
  int nsliceEff=nslice/nstride;
  fwd = fftw_plan_many_dft(1,&nsliceEff,nbin,
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_FORWARD,FFTW_MEASURE);
  fwd2 = fftw_plan_many_dft(1,&nsliceEff,1,
                            ptr2,0,1,nsliceEff,
                            ptr2,0,1,nsliceEff,
                            FFTW_FORWARD,FFTW_MEASURE);
  for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge; 
}

VIndEstimator::~VIndEstimator() {
  fftw_destroy_plan(fwd);
  fftw_destroy_plan(fwd2);
}

void VIndEstimator::initCalc(const int lnslice, const int firstSlice) {
  temp=0; temp2=0;
}

void VIndEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  int ibin=((int)(end[0]*dxinv+nbin))%nbin;
  int jbin=((int)(start[0]*dxinv+nbin))%nbin;
  if (ipart==0) {
    double sInd=coulAction->getAction(paths,islice);
    temp2(islice/nstride)+=sInd/tau;
  }
  if (ibin!=jbin) {
    int nstep = ((ibin-jbin+3*nbin/2)%nbin)-nbin/2;
    int idir = (nstep>0)?1:-1;
    if (idir>0) {
      for (int i=1; i<=nstep; i++) {
        temp((jbin+i)%nbin,islice/nstride)+=idir;
      }
    } else {
      for (int i=1; i<=-nstep; i++) {
        temp((ibin+i)%nbin,islice/nstride)+=idir;
      }
    }
  }
}

void VIndEstimator::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  blitz::Range allBin = blitz::Range::all();
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Reduce(&temp(0,0),&buff(0,0),
                                2*nslice/nstride*nbin,MPI::DOUBLE,MPI::SUM,0);
    temp(allBin,allSlice)=buff(allBin,allSlice); 
    mpi->getWorkerComm().Reduce(&temp2(0),&buff(0,0),
                                2*nslice/nstride,MPI::DOUBLE,MPI::SUM,0);
    temp2(allSlice)=temp2(allSlice); 
  }
#endif
  // Calculate autocorrelation function using FFT's.
  if (workerID==0) {
    fftw_execute(fwd);
    fftw_execute(fwd2);
    double betaInv=1./(tau*nslice);
    temp2*=tau;
    for (int ibin=0; ibin<nbin; ++ibin) {
      //temp(ibin,allSlice)*=conj(temp2(allSlice));
      for (int ifreq=0; ifreq<nfreq && 2*ifreq<nslice/2; ++ifreq) {
        value(ibin,ifreq) += real(temp(ibin,ifreq)*temp(ibin,ifreq)
                                  *conj(temp2(ifreq*2)))*betaInv;
      }
    }
    //value += imag(temp(allBin,blitz::Range(0,nfreq-1)))/(tau*nslice);
    norm+=1;
  }
}
