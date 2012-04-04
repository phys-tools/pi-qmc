// $Id: ExchangeMover.cc $
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
#include "ExchangeMover.h"
#include "Paths.h"
#include "Beads.h"
#include "util/RandomNumGenerator.h"
#include "util/SuperCell.h"
#include "SimulationInfo.h"
#include <cmath>
#include <sstream>
#include <string>

ExchangeMover::ExchangeMover(Paths& paths, const int& iFirstSlice,
		 const Vec dist, const MPIManager* mpi)
  : UniformMover(dist,mpi), paths(paths),
    iFirstSlice(iFirstSlice) {
}

ExchangeMover::~ExchangeMover() {
 }

double ExchangeMover::makeMove(VArray& displacement, const IArray& movingIndex) const {
  const Vec& bead1 = paths(movingIndex(0), iFirstSlice);
  const Vec& bead2 = paths(movingIndex(1), iFirstSlice);
  displacement(0) = bead2-bead1;
  displacement(1) = bead1-bead2;
#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    mpi->getWorkerComm().Bcast(displacement.data(),2*NDIM,MPI::DOUBLE,0);
 }
#endif
  return 0; 
}

