// $Id: MultiLevelSampler.cc 22 2009-5-18 Saad khairallah$
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
#include "stats/MPIManager.h"
#include "DoubleDisplaceMoveSampler.h"
#include "stats/AccRejEstimator.h"
#include "Action.h"
#include "DoubleAction.h"
#include "UniformMover.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include "Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include "SuperCell.h"
#include <sstream>
#include <string>

DoubleDisplaceMoveSampler::
DoubleDisplaceMoveSampler(int nmoving, int nrepeat,
			  Paths& paths, ParticleChooser& particleChooser, 
			  const UniformMover& mover, 
			  Action* action, DoubleAction* doubleAction, 
			  const MPIManager* mpi)
  : DisplaceMoveSampler(nmoving, nrepeat, paths, 
			particleChooser, mover, 
                        action, mpi),
    doubleAction(doubleAction),
    nsliceOver2(paths.getNSlice()/2) {
  iFirstSlice=(paths.getLowestOwnedSlice(true));
  iLastSlice=(paths.getHighestOwnedSlice(true));
  std::cout << "In DoubleDisplaceMove :: reassign iFirstSlice = " << iFirstSlice << std::endl;
  std::cout << "In DoubleDisplaceMove :: reassign iLastSlice  = " << iLastSlice << std::endl;
  std::cout << "In DoubleDisplaceMove :: nsliceOver2 = " << nsliceOver2 << std::endl;
}

DoubleDisplaceMoveSampler::~DoubleDisplaceMoveSampler() {}


bool DoubleDisplaceMoveSampler::tryMove() {
  int workerID = (mpi) ? mpi->getWorkerID() : 0;
  if (workerID==0) accRejEst->tryingMove(0);

  mover.makeMove(displacement,nmoving);
  if (nworker>1) {
    handleBoundary(iFirstSlice-1, iLastSlice+1, +1);
    handleBoundary(iFirstSlice-1+nsliceOver2, iLastSlice+1+nsliceOver2, +1);
  }
  // Evaluate the change in action.
  double deltaAction = (action==0) ?  0 
    : action->getActionDifference(paths,displacement,nmoving,movingIndex,
				  iFirstSlice,iLastSlice); 
  //Now do the other double section
  deltaAction += (action==0) ?  0 
    : action->getActionDifference(paths,displacement,nmoving,movingIndex,
				  iFirstSlice+nsliceOver2,iLastSlice+nsliceOver2);
  //Now consider nodal action
  deltaAction += (action==0) ?  0 
    : doubleAction->getActionDifference(paths,displacement,nmoving,movingIndex,
					iFirstSlice,iLastSlice);
  
#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    double netDeltaAction=0;
    mpi->getWorkerComm().Reduce(&deltaAction,&netDeltaAction,
                                1,MPI::DOUBLE,MPI::SUM,0);
    double acceptProb = exp(-netDeltaAction);
    bool acceptReject = RandomNumGenerator::getRand()>acceptProb;
    mpi->getWorkerComm().Bcast(&acceptReject,sizeof(bool),MPI::CHAR,0); 
    if (acceptReject) { 
      if (nworker>1) {
	handleBoundary(iFirstSlice-1, iLastSlice+1, -1);
	handleBoundary(iFirstSlice-1+nsliceOver2, iLastSlice+1+nsliceOver2, -1);
      }
      return false;
    }
  } else {
    double acceptProb=exp(-deltaAction); 
    if (RandomNumGenerator::getRand()>acceptProb)  { 
      if (nworker>1) {
	handleBoundary(iFirstSlice-1, iLastSlice+1, -1);
	handleBoundary(iFirstSlice-1+nsliceOver2, iLastSlice+1+nsliceOver2, -1);
      }
      return false;
    }
  }
#else
  double acceptProb=exp(-deltaAction);
    if (RandomNumGenerator::getRand()>acceptProb)  { 
      if (nworker>1) {
	handleBoundary(iFirstSlice-1, iLastSlice+1, -1);
	handleBoundary(iFirstSlice-1+nsliceOver2, iLastSlice+1+nsliceOver2, -1);
      }
      return false;
    }
#endif
  
  if (workerID==0) accRejEst->moveAccepted(0);
  
  // Move accepted.
  action->acceptLastMove();
  
  // Put moved beads in paths beads.
  for (int islice=iFirstSlice; islice<=iLastSlice; ++islice) {
    for (int imoving=0; imoving<nmoving; ++imoving) {
      Vec &bead(paths(movingIndex(imoving),islice));
      bead += displacement(imoving);
      bead = cell.pbc(bead);
  
      Vec &bead2(paths(movingIndex(imoving),islice+nsliceOver2));
      bead2 += displacement(imoving);
      bead2 = cell.pbc(bead2);
    }
  }
  
  return true;
}
