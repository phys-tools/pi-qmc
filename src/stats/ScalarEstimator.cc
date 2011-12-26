//$Id$
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
#include "ScalarEstimator.h"
#include "MPIManager.h"

ScalarEstimator::ScalarEstimator(const std::string& name)
  : Estimator(name,"","scalar"), scale(1.), shift(0.) {
}

ScalarEstimator::ScalarEstimator(const std::string &name,
  const std::string &typeString, const std::string &unitName,
  double scale, double shift)
  : Estimator(name,typeString,unitName), scale(scale), shift(shift) {
}

void ScalarEstimator::averageOverClones(const MPIManager* mpi) {
  int rank=0, size=1;
#ifdef ENABLE_MPI
  if (mpi->isCloneMain()) {
    rank = mpi->getCloneComm().Get_rank();
    size = mpi->getCloneComm().Get_size();
  }
#endif
  double v=calcValue(),value=v;
  reset();
#ifdef ENABLE_MPI
  if (mpi->isCloneMain()) {
    mpi->getCloneComm().Reduce(&v,&value,1,MPI::DOUBLE,MPI::SUM,0);
  }
#endif
  if (rank==0) setValue(value/size);
}
