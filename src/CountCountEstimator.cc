// $Id$
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
#include "CountCountEstimator.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <blitz/array.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Action.h"
#include "Paths.h"
#include "Distance.h"

CountCountEstimator::CountCountEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species *spec,
    const Vec &min, const Vec &max, const IVec &nbin, const IVecN &nbinN,
    const DistArray &dist, int nstride, MPIManager *mpi) 
  : BlitzArrayBlkdEst<2*NDIM+3>(name,nbinN,false),
    nslice(simInfo.getNSlice()), nfreq(nbinN[2*NDIM+2]), 
    nstride(nstride), maxc(nbinN[2*NDIM]), ntot(product(nbin)),
    min(min), deltaInv(nbin/(max-min)), nbin(nbin), nbinN(nbinN), dist(dist),
    cell(*simInfo.getSuperCell()), tau(simInfo.getTau()), temp(),
    count(nbin), ifirst(spec->ifirst), npart(spec->count), mpi(mpi) {
  scale=new VecN(1.);
  origin=new VecN(0.);
  for (int i=0; i<NDIM; ++i) {
    (*scale)[i] = (*scale)[i+NDIM] = (max[i]-min[i])/nbin[i];
    (*origin)[i] = (*origin)[i+NDIM] = min[i];
  }
  (*scale)[2*NDIM+2] = 2*3.141592653*simInfo.getTemperature();
  value=0.;
  blitz::TinyVector<int,NDIM+1> tempDim;
  for (int i=0; i<NDIM; ++i) tempDim[i]=nbin[i];
  tempDim[2*NDIM+2]=nslice/nstride;
  temp.resize(tempDim);
  // Set up new views of the arrays for convenience.
  count2 = new blitz::Array<int,1>(count.data(),
    blitz::shape(ntot), blitz::neverDeleteData); 
  temp2 = new blitz::Array<Complex,3>(temp.data(),
    blitz::shape(ntot,maxc,nslice/nstride), blitz::neverDeleteData); 
  value2 = new blitz::Array<float,5>(value.data(),
    blitz::shape(ntot,ntot,maxc,maxc,nfreq), blitz::neverDeleteData); 
  // Set up the FFT.
  fftw_complex *ptr = (fftw_complex*)temp.data();
  int nsliceEff=nslice/nstride;
  fwd = fftw_plan_many_dft(1,&nsliceEff,ntot*maxc,
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_FORWARD,FFTW_MEASURE);
}

CountCountEstimator::~CountCountEstimator() {
  for (int i=0; i<NDIM; ++i) delete dist[i];
  delete scale;
  delete origin;
  delete count2;
  delete temp2;
  delete value2;
  fftw_destroy_plan(fwd);
}

void CountCountEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0;
  count=0;
}


void CountCountEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec r=start;
    IVec ibin=0;
    for (int i=0; i<NDIM; ++i) {
      double d=(*dist[i])(r);
      ibin[i]=int(floor((d-min[i])*deltaInv[i]));
      if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
      if (i==NDIM-1) ++count(ibin);
    }
  }
  if (ipart==ifirst+npart-1) {
    for (int i=0; i<ntot; ++i) {
      int c = (*count2)(i);
      if (c<maxc) (*temp2)(i,c,islice/nstride) += 1.;
    }
    count=0;
  }
}


void CountCountEstimator::endCalc(const int nslice) {
  temp/=nstride;
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
    double betaInv=1./(tau*nslice);
    for (int i=0; i<ntot; ++i) 
      for (int j=0; j<ntot; ++j) 
        for (int ic=0; ic<maxc; ++ic) 
          for (int jc=0; jc<maxc; ++jc) 
            for (int ifreq=0; ifreq<nfreq; ++ifreq)
              (*value2)(i,j,ic,jc,ifreq) += betaInv
                *real((*temp2)(i,ic,ifreq)*conj((*temp2)(j,jc,ifreq)));
    norm+=1;
  }
}

void CountCountEstimator::reset() {}

void CountCountEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
