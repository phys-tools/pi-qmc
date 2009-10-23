// $Id$
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

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
#include "stats/MPIManager.h"
#include "RingLattice.h"
#include "Paths.h"
#include "Species.h"
#include <blitz/tinyvec-et.h>
#include "SuperCell.h"
#include "math.h"

RingLattice::RingLattice(Paths& paths, const double radius, const double angle0, MPIManager *mpi):paths(paths), radius(radius), angle0(angle0),ifirst(0),npart(paths.getNPart()), mpi(mpi) {
}

RingLattice::RingLattice(Paths& paths, const double radius, const double angle0, const Species& species, MPIManager *mpi):paths(paths), radius(radius), angle0(angle0), ifirst(species.ifirst), npart(species.count), mpi(mpi) {
}

void RingLattice::run() {
  int nslice = paths.getNSlice();
//  int nsites = npart;
//  IArray n(nsites);
//  n = npart / nsites;
//#ifdef ENABLE_MPI
//  if (mpi) mpi->getWorkerComm().Bcast(n.data(), nsites, MPI::INT, 0);
//#endif
  // Place the particles on the ring.
  int ifirstSlice = paths.getLowestSampleSlice(0,false);
  int ilastSlice = paths.getHighestStoredSlice(0,false);
  SuperCell cell=paths.getSuperCell();
  if (!mpi || mpi->isCloneMain()) {
    int ipart = ifirst;
    double angleDist = 2. * 3.141592653589793 / npart;
    for (int i=0; i<npart; ++i) {
      paths(ipart, ifirstSlice) = Vec(radius * cos(angleDist * i + angle0),
                                      radius * sin(angleDist * i + angle0));
      cell.pbc(paths(ipart, ifirstSlice));
      ++ipart;
    }
  }
  // Copy first slice to workers.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(&paths(ifirst,ifirstSlice), npart*NDIM,
                               MPI::DOUBLE,0);
  }
#endif
  // Copy to other slices.
  for (int islice=ifirstSlice+1; islice<ilastSlice; ++islice) {
    for (int ipart=ifirst; ipart<ifirst+npart; ++ipart) {
      paths(ipart,islice)=paths(ipart,ifirstSlice);
    }
  }
  if (paths.isDouble()) {
    for (int islice=ifirstSlice; islice<ilastSlice; ++islice) {
      for (int ipart=ifirst; ipart<ifirst+npart; ++ipart) {
        paths(ipart,islice+nslice/2)=paths(ipart,islice);
      }
    }
  }
}
