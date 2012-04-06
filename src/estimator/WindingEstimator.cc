// $Id$
/*  Copyright (C) 2009, 2011 John B. Shumway, Jr.

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
#include "util/SuperCell.h"
#include "Species.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include <fstream>
#include <blitz/tinyvec-et.h>

WindingEstimator::WindingEstimator(const SimulationInfo& simInfo,
  int nmax, const std::string &name, bool isChargeCoupled, const Species *species, MPIManager *mpi)
  : BlitzArrayBlkdEst<NDIM>(name,"histogram/winding", IVecN(2*nmax+1), true), 
    nmax(nmax), npart(simInfo.getNPart()), cell(*simInfo.getSuperCell()),
    charge(npart), isChargeCoupled(isChargeCoupled), ifirst(0), mpi(mpi) {
  value = 0.;
  norm = 0;
  charge = 1;
  if(species!=0) {
    ifirst = species->ifirst;
    npart = species->count;
  }

  if (isChargeCoupled) {
    for (int i=ifirst; i<npart; ++i) {
      charge(i)=int(simInfo.getPartSpecies(i).charge);
    }
  }
}

WindingEstimator::~WindingEstimator() {
}

void WindingEstimator::evaluate(const Paths &paths) {
  paths.sumOverLinks(*this);
}

void WindingEstimator::initCalc(const int nslice, const int firstSlice) {
  winding = 0.;
};

void WindingEstimator::handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec delta = start-end;
    cell.pbc(delta);
    winding += delta*charge(ipart);
  }
}

void WindingEstimator::endCalc(int nslice) {
  #ifdef ENABLE_MPI
  if (mpi) {
    Vec buffer;
    mpi->getWorkerComm().Reduce(&winding,&buffer,NDIM,MPI::DOUBLE,MPI::SUM,0);
    winding = buffer;
  }
  #endif
  IVec iwind = rint(winding*cell.b)+nmax;
  for (int idim=0; idim<NDIM; ++idim) {
    if (iwind[idim]<0 || iwind[idim]>2*nmax+1) break;
    if (idim==NDIM-1) ++value(iwind);
  }
  ++norm;
}

