// $Id: PermutationEstimator.cc,v 1.1 2008/01/08 20:46:40 jshumwa Exp $
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#include "PermutationEstimator.h"
#include "SimulationInfo.h"
#include "Paths.h"
#include "Permutation.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include <fstream>

PermutationEstimator::PermutationEstimator(const SimulationInfo& simInfo,
    MPIManager *mpi)
  : BlitzArrayBlkdEst<1>("permutations", IVecN(simInfo.getNPart()), true), 
    npart(simInfo.getNPart()), workFlags(npart),
    mpi(mpi) {
  value = 0.;
  norm = 0;
}

PermutationEstimator::~PermutationEstimator() {
}

void PermutationEstimator::evaluate(const Paths &paths) {
  const Permutation &perm(paths.getPermutation());
  if (perm[0]==0) {
    value(0)+=1.;
  } else {
    value(1)+=1.;
  }
  norm += 1.;
}
