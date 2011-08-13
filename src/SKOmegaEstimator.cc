// $Id: SKOmegaEstimator.cc 394 2011-08-05 02:10:09Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr.

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
#include "SKOmegaEstimator.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Action.h"
#include "Paths.h"

SKOmegaEstimator::SKOmegaEstimator(const SimulationInfo& simInfo,
    const std::string& name, const IVec &nbin, const IVecN &nbinN,
    int nstride, MPIManager *mpi) 
  : BlitzArrayBlkdEst<NDIM+3>(name,"dynamic-array/structure-factor",
                               nbinN,false),
    nslice(simInfo.getNSlice()), nfreq(nbinN[NDIM+2]), 
    nstride(nstride), npart(simInfo.getNPart()),
    nspec(simInfo.getNSpecies()),
    ntot(product(nbin)),
    min(-simInfo.getSuperCell()->a*0.5),
    deltaInv(-0.5*nbin/min), nbin(nbin), 
    cell(*simInfo.getSuperCell()), tau(simInfo.getTau()), temp(),
    speciesIndex(npart), mpi(mpi) {
  value=0.;
  for (int i=0; i<npart; ++i) {
    const Species *spec  = &simInfo.getPartSpecies(i);
    for (int j=0; j<nspec; ++j) {
      if (&simInfo.getSpecies(j) == spec) {
        speciesIndex(i) = j; break;
      }
    } 
  }
  int nsliceEff=nslice/nstride;
  blitz::TinyVector<int,NDIM+2> tempDim;
  tempDim[0] = nspec;
  for (int i=0; i<NDIM; ++i) tempDim[i+1] = nbin[i];
  tempDim[NDIM+1] = nsliceEff;
  temp.resize(tempDim);
  // Set up new views of the arrays for convenience.
  temp_ = new blitz::Array<Complex,3>(temp.data(),
    blitz::shape(nspec,ntot,nsliceEff), blitz::neverDeleteData); 
  value_ = new blitz::Array<float,4>(value.data(),
    blitz::shape(nspec,nspec,ntot,nfreq), blitz::neverDeleteData); 
  // Set up the FFT.
  fftw_complex *ptr = (fftw_complex*)temp.data();
  int n[NDIM+1];
  for (int i=0; i<NDIM; ++i) n[i]=nbin[i];
  n[NDIM]=nsliceEff; 
  fwd = fftw_plan_many_dft(NDIM+1,n,nspec,ptr,n,1,ntot*nsliceEff,
                                          ptr,n,1,ntot*nsliceEff,
                                          FFTW_FORWARD,FFTW_MEASURE);
}

SKOmegaEstimator::~SKOmegaEstimator() {
  delete temp_;
  delete value_;
  fftw_destroy_plan(fwd);
}

void SKOmegaEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0.;
}


void SKOmegaEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (islice%nstride==0) {
    int neffSlice = nslice/nstride;
    int isliceBin = (islice/nstride+neffSlice)%neffSlice;
    Vec r=start-min;
    blitz::TinyVector<int,NDIM+2> ibin=0;
    ibin[0] = speciesIndex(ipart);
    ibin[NDIM+1]=isliceBin;
    for (int i=0; i<NDIM; ++i) {
      ibin[i+1] = int(floor(r[i]*deltaInv[i]));
      if (ibin[i+1]<0 || ibin[i+1]>=nbin[i]) break;
      if (i==NDIM-1) temp(ibin) += 1.0;
    }
  }
}


void SKOmegaEstimator::endCalc(const int nslice) {
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
  ///Need code for multiple workers!
#ifdef ENABLE_MPI
//    if (mpi) {
//      if (workerID==0) {
//        mpi->getWorkerComm().Reduce(&temp(0,0,0),&temp(0,0,0),
//                                    product(nbin),MPI::DOUBLE,MPI::SUM,0);
//      } else {
//        mpi->getWorkerComm().Reduce(MPI::IN_PLACE,&temp(0,0,0),
//                                    product(nbin),MPI::DOUBLE,MPI::SUM,0);
//      }
//    }
#endif
  if (workerID==0) {
    // Calculate autocorrelation function using FFT for convolution.
    fftw_execute(fwd);
    double scale= tau*nstride*nstride/nslice;
    for (int ispec=0; ispec<nspec; ++ispec)
      for (int jspec=0; jspec<nspec; ++jspec)
        for (int k=0; k<ntot; ++k) 
          for (int kfreq=0; kfreq<nfreq; ++kfreq) 
          (*value_)(ispec,jspec,k,kfreq) += scale
            * real((*temp_)(ispec,k,kfreq)*conj((*temp_)(jspec,k,kfreq)));
    norm+=1;
  }
}

void SKOmegaEstimator::reset() {}

void SKOmegaEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
