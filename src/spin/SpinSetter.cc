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
#include "util/RandomNumGenerator.h"
#include "stats/MPIManager.h"
#include "SpinSetter.h"
#include "Paths.h"
#include "Species.h"
#include <cstdlib>
#include <blitz/tinyvec-et.h>
#include "util/SuperCell.h"

SpinSetter::SpinSetter(Paths& paths,MPIManager *mpi)
   : paths(paths), npart(paths.getNPart()), mpi(mpi) {
}

void SpinSetter::run() {
std::cout << "setting spin" << std::endl;
  int nslice=paths.getNSlice();
  int ifirstSlice=paths.getLowestOwnedSlice(false);
  int ilastSlice=paths.getHighestOwnedSlice(false);
  SuperCell cell=paths.getSuperCell();
  if (!mpi || mpi->isCloneMain()) {
    for (int ipart=0; ipart<npart; ++ipart) {
      (*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,ifirstSlice,1)))
      =SVec(-0.485332,-0.433568,0.476574,0.375401);
    }
  }
  // Copy first slice to workers.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(paths.getAuxBead(0,ifirstSlice,1), npart*4,
                               MPI::DOUBLE,0);
  } 
#endif
  // Copy to other slices.
  for (int islice=ifirstSlice+1; islice<ilastSlice; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      (*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,islice,1)))
        =(*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,ifirstSlice,1)));
    } 
  }
  if (paths.isDouble()) {
  for (int islice=ifirstSlice; islice<ilastSlice; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      (*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,islice+nslice/2,1)))
        =(*reinterpret_cast<SVec*>(paths.getAuxBead(ipart,islice,1)));
    } 
  }
  }
}
