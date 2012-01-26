// $Id: DensityCurrentEstimator.cc 394 2011-08-05 02:10:09Z john.shumwayjr $
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#include "DensityCurrentEstimator.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Action.h"
#include "Paths.h"
#include "Distance.h"

DensityCurrentEstimator::DensityCurrentEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Vec &min, const Vec &max, const IVec &nbin, 
    const IVecN &nbinN, const DistArray &dist, int nstride, MPIManager *mpi) 
  : BlitzArrayBlkdEst<NDIM+1>(name,"dynamic-array/density-current",
                               nbinN,true),
//    nslice(simInfo.getNSlice()), tempn(),
    nsliceEff(simInfo.getNSlice()/nstride), nfreq(nbinN[NDIM]), 
    nstride(nstride), min(min), deltaInv(nbin/(max-min)), nbin(nbin), 
    tempj(1,simInfo.getNSlice()/nstride), dist(dist), tau(simInfo.getTau()), 
    ax(simInfo.getSuperCell()->a[0]/2.), ntot(product(nbin)),  
    npart(simInfo.getNPart()), q(npart), mpi(mpi) {
  DensityCurrentEstimator::IVecN tempDim;
  for (int i=0; i<NDIM; ++i) tempDim[i]=nbin[i];
  tempDim[NDIM]=nsliceEff;
  tempn.resize(tempDim);
  // Set up new views of the arrays for convenience.
  tempn_ = new DensityCurrentEstimator::CArray2(tempn.data(),
    blitz::shape(ntot,nsliceEff), blitz::neverDeleteData); 
  value_ = new DensityCurrentEstimator::FArray2(value.data(),
    blitz::shape(ntot,nfreq), blitz::neverDeleteData); 
  // Set up the FFT.
  fftw_complex *ptr = (fftw_complex*)tempn.data();
  fwdn = fftw_plan_many_dft(1,&nsliceEff,ntot,
                            ptr,0,1,nsliceEff,
                            ptr,0,1,nsliceEff,
                            FFTW_FORWARD,FFTW_MEASURE);
  ptr = (fftw_complex*)tempj.data();
  int tot = 1;
  fwdj = fftw_plan_many_dft(1,&nsliceEff,tot,
                            ptr,0,1,nsliceEff,
                            ptr,0,1,nsliceEff,
                            FFTW_FORWARD,FFTW_MEASURE);
  for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge; 
}

DensityCurrentEstimator::~DensityCurrentEstimator() {
  for (int i=0; i<NDIM; ++i) delete dist[i];
  delete tempn_;
  delete value_;
  fftw_destroy_plan(fwdn);
  fftw_destroy_plan(fwdj);
}

void DensityCurrentEstimator::initCalc(const int nslice,
    const int firstSlice) {
  tempn=0.;
  tempj=0.;
}


void DensityCurrentEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  // Calculate density distribution.
  int isliceBin = (islice/nstride+nsliceEff)%nsliceEff;
  Vec r=start;
  blitz::TinyVector<int,NDIM+1> ibin=0;
  ibin[NDIM]=isliceBin;
  for (int i=0; i<NDIM; ++i) {
    double d=(*dist[i])(r);
    ibin[i]=int(floor((d-min[i])*deltaInv[i]));
    if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
    if (i==NDIM-1) tempn(ibin) += 1.0;
  }
  // Calculate current at x = 0, assuming no link is longer than a[0]/2.
  if (fabs(start[0]-end[0]) < ax && start[0]*end[0] < 0.) {
    if (start[0] > end[0]) tempj(0,isliceBin) -= q(ipart);
    else tempj(0,isliceBin) += q(ipart);
  }
}


void DensityCurrentEstimator::endCalc(const int nslice) {
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
    fftw_execute(fwdn);
    fftw_execute(fwdj);
//    double scale= tau/nsliceEff;
    for (int i=0; i<ntot; ++i) 
      for (int ifreq=0; ifreq<nfreq; ++ifreq)
	(*value_)(i,ifreq) += //scale *
	  real((*tempn_)(i,ifreq)*conj(tempj(0,ifreq)));
    norm+=1;
  }
}

