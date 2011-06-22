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
#include "DensDensEstimator.h"
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

DensDensEstimator::DensDensEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species *spec1, const Species *spec2,
    const Vec &min, const Vec &max, const IVec &nbin, const IVecN &nbinN,
    const DistArray &dist, int nstride, MPIManager *mpi) 
  : BlitzArrayBlkdEst<2*NDIM+1>(name,"dynamic-array/density-density",
                               nbinN,false),
    nslice(simInfo.getNSlice()), nfreq(nbinN[2*NDIM]), 
    nstride(nstride), ntot(product(nbin)),
    min(min), deltaInv(nbin/(max-min)), nbin(nbin), nbinN(nbinN), dist(dist),
    cell(*simInfo.getSuperCell()), tau(simInfo.getTau()), temp1(), temp2(),
    ifirst1(spec1->ifirst), npart1(spec1->count), 
    ifirst2(spec2->ifirst), npart2(spec2->count), mpi(mpi) {
  scale=new VecN(1.);
  origin=new VecN(0.);
  for (int i=0; i<NDIM; ++i) {
    (*scale)[i] = (*scale)[i+NDIM] = (max[i]-min[i])/nbin[i];
    (*origin)[i] = (*origin)[i+NDIM] = min[i];
  }
  (*scale)[2*NDIM] = 2*3.141592653*simInfo.getTemperature();
  value=0.;
  blitz::TinyVector<int,NDIM+1> tempDim;
  for (int i=0; i<NDIM; ++i) tempDim[i]=nbin[i];
  tempDim[NDIM]=nslice/nstride;
  temp1.resize(tempDim);
  temp2.resize(tempDim);
  // Set up new views of the arrays for convenience.
  temp1_ = new blitz::Array<Complex,2>(temp1.data(),
    blitz::shape(ntot,nslice/nstride), blitz::neverDeleteData); 
  temp2_ = new blitz::Array<Complex,2>(temp2.data(),
    blitz::shape(ntot,nslice/nstride), blitz::neverDeleteData); 
  value_ = new blitz::Array<float,3>(value.data(),
    blitz::shape(ntot,ntot,nfreq), blitz::neverDeleteData); 
  // Set up the FFT.
  fftw_complex *ptr = (fftw_complex*)temp1.data();
  int nsliceEff=nslice/nstride;
  fwd1 = fftw_plan_many_dft(1,&nsliceEff,ntot,
                            ptr,0,1,nsliceEff,
                            ptr,0,1,nsliceEff,
                            FFTW_FORWARD,FFTW_MEASURE);
  ptr = (fftw_complex*)temp2.data();
  fwd2 = fftw_plan_many_dft(1,&nsliceEff,ntot,
                            ptr,0,1,nsliceEff,
                            ptr,0,1,nsliceEff,
                            FFTW_FORWARD,FFTW_MEASURE);
}

DensDensEstimator::~DensDensEstimator() {
  for (int i=0; i<NDIM; ++i) delete dist[i];
  delete scale;
  delete origin;
  delete temp1_;
  delete temp2_;
  delete value_;
  fftw_destroy_plan(fwd1);
  fftw_destroy_plan(fwd2);
}

void DensDensEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp1=0;
  temp2=0;
}


void DensDensEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (islice%nstride==0 && ipart>=ifirst1 && ipart<ifirst1+npart1) {
    Vec r=start;
    blitz::TinyVector<int,NDIM+1> ibin=0;
    ibin[NDIM]=islice/nstride;
    for (int i=0; i<NDIM; ++i) {
      double d=(*dist[i])(r);
      ibin[i]=int(floor((d-min[i])*deltaInv[i]));
      if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
      if (i==NDIM-1) temp1(ibin) += 1.0;
    }
  }
  if (ipart>=ifirst2 && ipart<ifirst2+npart2) {
    Vec r=start;
    blitz::TinyVector<int,NDIM+1> ibin=0;
    ibin[NDIM]=islice/nstride;
    for (int i=0; i<NDIM; ++i) {
      double d=(*dist[i])(r);
      ibin[i]=int(floor((d-min[i])*deltaInv[i]));
      if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
      if (i==NDIM-1) temp2(ibin) += 1.0;
    }
  }
}


void DensDensEstimator::endCalc(const int nslice) {
  temp1/=nstride;
  temp2/=nstride;
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
    fftw_execute(fwd1);
    fftw_execute(fwd2);
    double betaInv=1./(tau*nslice);
    for (int i=0; i<ntot; ++i) 
      for (int j=0; j<ntot; ++j) 
        for (int ifreq=0; ifreq<nfreq; ++ifreq)
          (*value_)(i,j,ifreq) += betaInv
            *real((*temp1_)(i,ifreq)*conj((*temp2_)(j,ifreq)));
    norm+=1;
  }
}

void DensDensEstimator::reset() {}

void DensDensEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
