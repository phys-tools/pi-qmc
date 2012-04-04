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
#include "DynamicPCFEstimator.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "action/Action.h"
#include "Paths.h"
#include "util/PairDistance.h"

DynamicPCFEstimator::DynamicPCFEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species *spec1, const Species *spec2,
    double min, double max, int nbin, const int nfreq, int nstride,
    const PairDistance *dist, MPIManager *mpi) 
  : BlitzArrayBlkdEst<3>(name,"dynamic-array/pair-correlation",
                         IVecN(nbin,nbin,nfreq),false),
    nslice(simInfo.getNSlice()), nfreq(nfreq), nstride(nstride), nbin(nbin),
    min(min), deltaInv(nbin/(max-min)), dist(dist),
    cell(*simInfo.getSuperCell()), tau(simInfo.getTau()), temp(),
    ifirst(spec1->ifirst), jfirst(spec2->ifirst),
    nipart(spec1->count), njpart(spec2->count), mpi(mpi) {
  scale=new VecN(1.);
  origin=new VecN(0.);
  (*scale)[0] = (*scale)[1] = (max-min)/nbin;
  (*origin)[0] = (*origin)[1] = min;
  (*scale)[2] = 2*3.141592653*simInfo.getTemperature();
  value=0.;
  blitz::TinyVector<int,2> tempDim;
  tempDim[0]=nbin;
  tempDim[1]=nslice/nstride;
  temp.resize(tempDim);
  // Set up the FFT.
  fftw_complex *ptr = (fftw_complex*)temp.data();
  int nsliceEff=nslice/nstride;
  fwd = fftw_plan_many_dft(1,&nsliceEff,nbin,
                           ptr,0,1,nsliceEff,
                           ptr,0,1,nsliceEff,
                           FFTW_FORWARD,FFTW_MEASURE);
}

DynamicPCFEstimator::~DynamicPCFEstimator() {
  delete dist;
  delete scale;
  delete origin;
  fftw_destroy_plan(fwd);
}

void DynamicPCFEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0;
}


void DynamicPCFEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (ipart>=ifirst && ipart<ifirst+nipart) {
    Vec r1=start;
    for (int jpart=jfirst; jpart<(jfirst+njpart); ++jpart) {
      if (ipart!=jpart) {
        Vec r2=paths(jpart,islice);
        int ibin = 0;
        double d = (*dist)(r1,r2,cell);
        ibin = int(floor((d-min)*deltaInv));
        if (ibin>=0 && ibin<nbin) temp(ibin,islice/nstride) += 1.0;
      }
    }
  }
}


void DynamicPCFEstimator::endCalc(const int nslice) {
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
    for (int i=0; i<nbin; ++i) 
      for (int j=0; j<nbin; ++j) 
        for (int ifreq=0; ifreq<nfreq; ++ifreq)
          value(i,j,ifreq) += betaInv
            *real(temp(i,ifreq)*conj(temp(j,ifreq)));
    norm+=1;
  }
}

void DynamicPCFEstimator::reset() {}

void DynamicPCFEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
