// $Id: CollectiveMover.cc $
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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "stats/MPIManager.h"
#include "CollectiveMover.h"
#include "Paths.h"
#include "Beads.h"
#include "RandomNumGenerator.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include <cmath>

#include <sstream>
#include <string>

CollectiveMover::CollectiveMover(Paths& paths, const int& iFirstSlice,
                 const Vec kvec, const Vec amp, const MPIManager* mpi)
  : UniformMover(amp,mpi), paths(paths), amp(amp), kvec(kvec),
    iFirstSlice(iFirstSlice) {
  std::cout<<"CollectiveMover: Amplitude = "<<amp<<", k = "<<kvec<<std::endl;
}

CollectiveMover::~CollectiveMover() {
 }

double CollectiveMover::makeMove(VArray& displacement, const IArray& movingIndex) const {
  int npart = displacement.size();
  const Vec &referBead = paths(movingIndex(0), iFirstSlice);
//#ifdef ENABLE_MPI
//  if (mpi && (mpi->getNWorker())>1) {
//    mpi->getWorkerComm().Bcast(&(referBead[0]),NDIM,MPI::DOUBLE,0);
// }
//#endif
  for (int ipart=0;ipart<npart;++ipart) {
    displacement(ipart) = amp * 
           cos(dot(kvec, paths(ipart,iFirstSlice) - referBead));
  }
#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    mpi->getWorkerComm().Bcast(displacement.data(),npart*NDIM,MPI::DOUBLE,0);
 }
#endif
  return 0; 
}

