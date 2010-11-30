// $Id$
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
#include "ModelSampler.h"
#include "stats/AccRejEstimator.h"
#include "stats/MPIManager.h"
#include "Action.h"
#include "ActionChoice.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include <sstream>
#include <string>

ModelSampler::ModelSampler(Paths& paths, Action* action, 
  ActionChoiceBase* actionChoice, const MPIManager* mpi)
  : paths(paths), action(action), actionChoice(actionChoice),
    nmodel(actionChoice->getModelCount()), accRejEst(0),  mpi(mpi) {
}

ModelSampler::~ModelSampler() {
}

void ModelSampler::run() {
  tryMove();

}


bool ModelSampler::tryMove() {

  accRejEst->tryingMove(0);

  int imodel = actionChoice->getModelState(); 

  // We select the next model so that the end points always
  // try a valid model. This is efficient for two models, but
  // efficiency is around 60% for several models.
  // Probably want to revisit this later.
  int jmodel = imodel+2*((int)floor(RandomNumGenerator::getRand()
      +(nmodel-1.-imodel)/(nmodel-1.)))-1;
  if (jmodel<0 or jmodel>=nmodel) return false;

#ifdef ENABLE_MPI
  int nworker = (mpi) ? mpi->getNWorker():1;
  if (mpi && nworker>1) {
    mpi->getWorkerComm().Bcast(&jmodel,1,MPI::INT,0);
  }
#endif 

  double tranProb = (jmodel>imodel)
                    ? (nmodel-1.-imodel)/jmodel :imodel/(nmodel-1.-jmodel);

  // Evaluate the change in action.
  double deltaAction = actionChoice->getActionDifference(paths,jmodel);
#ifdef ENABLE_MPI
  double totalDeltaAction = 0;
  mpi->getWorkerComm().Reduce(&deltaAction,&totalDeltaAction,
                              1,MPI::DOUBLE,MPI::SUM,0);
  deltaAction = totalDeltaAction;
#endif 

  double acceptProb=exp(-deltaAction)/tranProb;

  bool reject = RandomNumGenerator::getRand()>acceptProb;
#ifdef ENABLE_MPI
    if (nworker > 1) {
      mpi->getWorkerComm().Bcast(&reject, sizeof(bool), MPI::CHAR, 0); 
    }
#endif 
  if (reject) return false;

  actionChoice->setModelState(jmodel);
  paths.setModelState(jmodel);
  accRejEst->moveAccepted(0);
  return true;
}

AccRejEstimator* 
ModelSampler::getAccRejEstimator(const std::string& name) {
  return accRejEst=new AccRejEstimator(name.c_str(),1);
}
