// $Id: WindingEstimator.cc 22 2009-03-06 13:52:07Z john.shumwayjr $
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
#include "WindingEstimator.h"
#include "SimulationInfo.h"
#include "Paths.h"
#include "SuperCell.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include <fstream>
#include <blitz/tinyvec-et.h>

WindingEstimator::WindingEstimator(const SimulationInfo& simInfo,
  int nmax, MPIManager *mpi)
  : BlitzArrayBlkdEst<4>("winding",
                         IVecN(2*nmax+1,2*nmax+1,2*nmax+1,2*nmax+1), true), 
    mpi(mpi), nmax(nmax), npart(simInfo.getNPart()),
    cell(*simInfo.getSuperCell()) {
  value = 0.;
  norm = 0;
}

WindingEstimator::~WindingEstimator() {
}

void WindingEstimator::evaluate(const Paths &paths) {
  paths.sumOverLinks(*this);
}

void WindingEstimator::initCalc(const int nslice, const int firstSlice) {
  winding1 = 0.;
  winding2 = 0.;
};

void WindingEstimator::handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&) {
  Vec delta = start-end;
  cell.pbc(delta);
  if (ipart<npart/2) {
    winding1 += delta;
  } else {
    winding2 += delta;
  }
}

void WindingEstimator::endCalc(int nslice) {
  #ifdef ENABLE_MPI
  if (mpi) {
    Vec buffer1, buffer2;
    mpi->getWorkerComm().Reduce(&winding1,&buffer1,NDIM,MPI::DOUBLE,MPI::SUM,0);
    mpi->getWorkerComm().Reduce(&winding2,&buffer2,NDIM,MPI::DOUBLE,MPI::SUM,0);
    winding1 = buffer1; winding2 = buffer2;
  }
  #endif
  IVec iwind1 = rint(winding1*cell.b);
  IVec iwind2 = rint(winding2*cell.b);
#if NDIM>1
  if (iwind1[0]>=-nmax && iwind1[0]<=nmax &&
      iwind1[1]>=-nmax && iwind1[1]<=nmax &&
      iwind2[0]>=-nmax && iwind2[0]<=nmax &&
      iwind2[1]>=-nmax && iwind2[1]<=nmax)
    ++value(iwind1[0]+nmax,iwind1[1]+nmax,iwind2[0]+nmax,iwind2[1]+nmax);
#endif
  ++norm;
}

