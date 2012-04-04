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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "WorkerShifter.h"
#include "Paths.h"
#include "util/RandomNumGenerator.h"
#include "stats/MPIManager.h"

WorkerShifter::WorkerShifter(const int maxShift, Paths& paths, MPIManager* mpi)
  : CompositeAlgorithm(1), maxShift(maxShift), paths(paths), mpi(mpi) {
}

WorkerShifter::~WorkerShifter() {}

void WorkerShifter::run() {
  int ishift=(int)((maxShift-1)*RandomNumGenerator::getRand()*(1-1e-8))+1;
  if (maxShift==0) ishift=0;
#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(&ishift,1,MPI::INT,0);
#endif
  paths.shift(ishift);
  CompositeAlgorithm::run();
}
