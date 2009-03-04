// $Id: DensityEstimator.cc 12 2009-02-07 23:32:51Z john.shumwayjr $
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
#include "DensityEstimator.h"
#include "stats/BlitzArrayBlkdEst.h"
#include "stats/MPIManager.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <blitz/array.h>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Action.h"
#include "Paths.h"

DensityEstimator::DensityEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species &spec, int N,
    const Vec &min, const Vec &max, const IVec &nbin,
    const DistArray &dist, MPIManager *mpi) 
  : ArrayBlockedEstimator(name,true),
    N(N), min(min), deltaInv(nbin/(max-min)), nbin(nbin), dist(dist),
    cell(*simInfo.getSuperCell()), 
    ifirst(spec.ifirst), npart(spec.count), mpi(mpi) {
  //BlitzArrayBlkdEst<N>::norm=0;
}


DensityEstimator::~DensityEstimator() {
  for (int i=0; i<N; ++i) delete dist[i];
}

void DensityEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0;
}


void DensityEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec r=start;
    //IVecN ibin=0;
    for (int i=0; i<N; ++i) {
      //double d=(*dist[i])(r,cell);
      //ibin[i]=int((d-min[i])*deltaInv[i]);
      //if (d<min[i] || ibin[i]>nbin[i]-1) break;
      //if (i==N-1) ++temp(ibin);
    }
  }
}


void DensityEstimator::endCalc(const int nslice) {
  temp/=nslice;
  //BlitzArrayBlkdEst<N>::norm+=1;
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
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
  ///Need code for multiple workers!
  if (workerID==0) {
    //BlitzArrayBlkdEst<N>::value+=temp;
  }
}

void DensityEstimator::reset() {}

void DensityEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
