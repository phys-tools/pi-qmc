//$Id: SpinChoiceFixedNodeAction.cc 383 2011-04-20 16:56:02Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr. and Jianheng Liu

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
#include "SpinChoiceFixedNodeAction.h"
#include "SpinModelState.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>

SpinChoiceFixedNodeAction::SpinChoiceFixedNodeAction(
  const SimulationInfo &simInfo, int initial,
  const Species &species, NodeModel *nodeModel, bool withNodalAction,
  bool useDistDerivative, int maxlevel, bool useManyBodyDistance,
  const MPIManager* mpi) 
  : FixedNodeAction(simInfo,species,nodeModel,withNodalAction,
      useDistDerivative,maxlevel,useManyBodyDistance), mpi(mpi) {
  std::cout << "npart for spin flip is " << nSpeciesPart << std::endl;
  spinModelState = new SpinModelState(nSpeciesPart,initial);
  modelState = spinModelState;
  nodeModel->setSpinModelState(spinModelState);
}

SpinChoiceFixedNodeAction::~SpinChoiceFixedNodeAction() {
  delete spinModelState;
}

void SpinChoiceFixedNodeAction::initCalc(const int nslice, 
    const int firstSlice) {
  totalAction = 0.;
}

double SpinChoiceFixedNodeAction::getActionDifference(const Paths &paths,
    int ipart) {
  paths.sumOverLinks(*this);
  double oldAction = totalAction;
  spinModelState->flipSpin(ipart);
#ifdef ENABLE_MPI
  if (mpi) {
//      mpi->getWorkerComm().Bcast(spinModelState->getModelState().data(),
//                                 spinModelState->getModelCount(),MPI::INT,0);
    spinModelState->broadcastToMPIWorkers(mpi);
  }
#endif
  paths.sumOverLinks(*this);
  double newAction = totalAction;
  spinModelState->flipSpin(ipart);
#ifdef ENABLE_MPI
  if (mpi) {
//    mpi->getWorkerComm().Bcast(spinModelState->getModelState().data(),
//                               spinModelState->getModelCount(),MPI::INT,0);
    spinModelState->broadcastToMPIWorkers(mpi);
  }       
#endif

//  std::cout << "newAction, oldAction " << newAction << ", " << oldAction << std::endl; 

  return newAction-oldAction;
}

void SpinChoiceFixedNodeAction::handleLink(const LinkSummable::Vec &start,
            const LinkSummable::Vec &end, int ipart, int islice, 
	    const Paths &paths) {
  double u=0., utau=0., ulambda=0;
  FixedNodeAction::Vec fm=0., fp=0.;
  this->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
  totalAction += u;
}

