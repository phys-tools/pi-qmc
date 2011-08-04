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
#include "DensityEstimator.h"
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

DensityEstimator::DensityEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species *spec,
    const Vec &min, const Vec &max, const IVec &nbin,
    const DistArray &dist, MPIManager *mpi) 
  : BlitzArrayBlkdEst<NDIM>(name,"array/density",nbin,false),
    min(min), deltaInv(nbin/(max-min)), nbin(nbin), dist(dist),
    cell(*simInfo.getSuperCell()),temp(nbin),
#ifdef ENABLE_MPI
    mpiBuffer(nbin),
#endif 
    ifirst(spec->ifirst), npart(spec->count), mpi(mpi) {
  scale=new Vec((max-min)/nbin);
  origin=new Vec(min);
  value=0.;
}


DensityEstimator::~DensityEstimator() {
  for (int i=0; i<NDIM; ++i) delete dist[i];
  delete scale;
  delete origin;
}

void DensityEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0.;
}


void DensityEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec r=start;
    IVec ibin=0;
    for (int i=0; i<NDIM; ++i) {
      double d=(*dist[i])(r);
      ibin[i]=int(floor((d-min[i])*deltaInv[i]));
      if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
      if (i==NDIM-1) ++temp(ibin);
    }
  }
}


void DensityEstimator::endCalc(const int lnslice) {
  int nslice = lnslice;
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    int ibuffer;
    mpi->getWorkerComm().Reduce(temp.data(),mpiBuffer.data(),
                                product(nbin),MPI::FLOAT,MPI::SUM,0);
    mpi->getWorkerComm().Reduce(&lnslice,&ibuffer,1,MPI::INT,MPI::SUM,0);
    temp = mpiBuffer;
    nslice = ibuffer;
  }
#endif
  temp /= nslice;
  if (workerID==0) {
    BlitzArrayBlkdEst<NDIM>::value+=temp;
    norm+=1.;
  }
}

void DensityEstimator::reset() {}

void DensityEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
