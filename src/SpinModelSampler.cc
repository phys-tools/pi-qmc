// $Id: SpinModelSampler.cc 404 2011-09-19 21:45:15Z john.shumwayjr $
/*  Copyright (C) 2010 John B. Shumway, Jr.

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
#include "SpinModelSampler.h"
#include "stats/AccRejEstimator.h"
#include "stats/MPIManager.h"
#include "Action.h"
#include "ActionChoice.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include "Permutation.h"
#include "SpinModelState.h"
#include <sstream>
#include <string>

SpinModelSampler::SpinModelSampler(Paths& paths, Action* action, 
  ActionChoiceBase* actionChoice, const MPIManager* mpi)
  : paths(paths), action(action), actionChoice(actionChoice),
    modelState(dynamic_cast<SpinModelState&>
                 (actionChoice->getModelState())),
    nmodel(modelState.getModelCount()), accRejEst(0),  mpi(mpi) 
#ifdef ENABLE_MPI
    ,npart(paths.getNPart()), nworker((mpi)?mpi->getNWorker():1),
    iworkerPerm((mpi)?mpi->getNWorker():1,paths.getNPart())
#endif
{
}

SpinModelSampler::~SpinModelSampler() {
}

void SpinModelSampler::run() {
  tryMove();

}


bool SpinModelSampler::tryMove() {
  int workerID = (mpi) ? mpi->getWorkerID() : 0;

  if (workerID==0) accRejEst->tryingMove(0);

  int ipart = nmodel;
  do {
    ipart  = RandomNumGenerator::getRand() * (nmodel-1);
//std::cout << ipart << ", " << nmodel << std::endl;
  } while (! (ipart < nmodel-1));

  // Check if ipart is part of a permutation, reject if yes.
  Permutation permutation = paths.getGlobalPermutation();
  if (ipart != permutation[ipart]) return false;

  // Evaluate the change in action.
  double deltaAction = actionChoice->getActionDifference(paths,ipart);

#ifdef ENABLE_MPI
  double totalDeltaAction = 0;
  mpi->getWorkerComm().Reduce(&deltaAction,&totalDeltaAction,
                              1,MPI::DOUBLE,MPI::SUM,0);
  deltaAction = totalDeltaAction;
#endif 

  double acceptProb=exp(-deltaAction);

  bool reject = RandomNumGenerator::getRand()>acceptProb;
#ifdef ENABLE_MPI
    if (nworker > 1) {
      mpi->getWorkerComm().Bcast(&reject, sizeof(bool), MPI::CHAR, 0); 
    }
#endif 
  if (reject) return false;

  modelState.flipSpin(ipart);

  if (workerID==0) accRejEst->moveAccepted(0);
  return true;
}

AccRejEstimator* 
SpinModelSampler::getAccRejEstimator(const std::string& name) {
  return accRejEst=new AccRejEstimator(name.c_str(),1);
}
