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
#include "NodeModel.h"
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

double SpinChoiceFixedNodeAction::getActionChoiceDifference(const Paths &paths,
    int ipart) {
  int totNSlice = paths.getNSlice();
  double newAction = 0., oldAction = 0.;
  int iFirstSlice = paths.getLowestOwnedSlice(false);
  int nSlice = paths.getHighestOwnedSlice(false);
  spinModelState->flipSpin(ipart);
#ifdef ENABLE_MPI
  if (mpi) {
    spinModelState->broadcastToMPIWorkers(mpi);
  }
#endif
  for (int islice=iFirstSlice; islice<=nSlice; ++islice) {
    int jslice = (islice+totNSlice/2)%totNSlice;
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r1(i)=paths(i,islice,-1);
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r2(i)=paths(i,jslice,-1);
    NodeModel::DetWithFlag detm = 
      FixedNodeAction::nodeModel->evaluate(FixedNodeAction::r1,
	                                   FixedNodeAction::r2, 0, false);
    FixedNodeAction::nodeModel->evaluateDistance(FixedNodeAction::r1,
	                                         FixedNodeAction::r2,
						 0, FixedNodeAction::dim1,
						 FixedNodeAction::dim2);
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r1(i)=paths(i,islice);
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r2(i)=paths(i,jslice);
    NodeModel::DetWithFlag det = 
      FixedNodeAction::nodeModel->evaluate(FixedNodeAction::r1,
                                           FixedNodeAction::r2, 0, false);
    FixedNodeAction::nodeModel->evaluateDistance(FixedNodeAction::r1,
	                                         FixedNodeAction::r2,0,
						 FixedNodeAction::di1,
						 FixedNodeAction::di2);
    // check if nodes are crossed.
    if (detm.err || detm.err || detm.det*det.det < 0) 
    {newAction += 1e200;continue;}

    // Calculate action.
    if (FixedNodeAction::useManyBodyDistance) {
      double d12=0, d1m2=0;
      for (int i=0; i<npart; ++i) {
	d12 += 1./(FixedNodeAction::di1(i)*FixedNodeAction::di1(i));
	d1m2 += 1./(FixedNodeAction::dim1(i)*FixedNodeAction::dim1(i));
      }
      d12 = 1./sqrt(d12);
      dim2 = 1./sqrt(d1m2);
      double xim1 = d12*d1m2;
      double exim1 = exp(-xim1);
      newAction -= log((1-exim1));
    } else {
      for (int jpart=0; jpart<npart; ++jpart) {
	double xim1 = FixedNodeAction::dim1(jpart)*FixedNodeAction::di1(jpart);
	newAction += -log(1-exp(-xim1));
      }
    }
  }
  // Then calculate the old action.
  spinModelState->flipSpin(ipart);
#ifdef ENABLE_MPI
  if (mpi) {
    spinModelState->broadcastToMPIWorkers(mpi);
  }
#endif
  for (int islice=iFirstSlice; islice<=nSlice; ++islice) {
    int jslice = (islice+totNSlice/2)%totNSlice;
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r1(i)=paths(i,islice,-1);
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r2(i)=paths(i,jslice,-1);
    NodeModel::DetWithFlag detm = 
      FixedNodeAction::nodeModel->evaluate(FixedNodeAction::r1,
	                                   FixedNodeAction::r2, 0, false);
    FixedNodeAction::nodeModel->evaluateDistance(FixedNodeAction::r1,
	                                         FixedNodeAction::r2,
						 0, FixedNodeAction::dim1,
						 FixedNodeAction::dim2);
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r1(i)=paths(i,islice);
    for (int i=0; i<npart; ++i) 
      FixedNodeAction::r2(i)=paths(i,jslice);
    NodeModel::DetWithFlag det = 
      FixedNodeAction::nodeModel->evaluate(FixedNodeAction::r1,
                                           FixedNodeAction::r2, 0, false);
    FixedNodeAction::nodeModel->evaluateDistance(FixedNodeAction::r1,
	                                         FixedNodeAction::r2,0,
						 FixedNodeAction::di1,
						 FixedNodeAction::di2);
    // check if nodes are crossed.
    if (detm.err || detm.err || detm.det*det.det < 0) 
    {oldAction += 1e200;continue;}

    // Calculate action.
    if (FixedNodeAction::useManyBodyDistance) {
      double d12=0, d1m2=0;
      for (int i=0; i<npart; ++i) {
	d12 += 1./(FixedNodeAction::di1(i)*FixedNodeAction::di1(i));
	d1m2 += 1./(FixedNodeAction::dim1(i)*FixedNodeAction::dim1(i));
      }
      d12 = 1./sqrt(d12);
      dim2 = 1./sqrt(d1m2);
      double xim1 = d12*d1m2;
      double exim1 = exp(-xim1);
      oldAction -= log((1-exim1));
    } else {
      for (int jpart=0; jpart<npart; ++jpart) {
	double xim1 = FixedNodeAction::dim1(jpart)*FixedNodeAction::di1(jpart);
	oldAction += -log(1-exp(-xim1));
      }
    }
  }
  return newAction - oldAction;
}

// void SpinChoiceFixedNodeAction::initCalc(const int nslice, 
//     const int firstSlice) {
//   totalAction = 0.;
// }
// 
// double SpinChoiceFixedNodeAction::getActionDifference(const Paths &paths,
//     int ipart) {
//   paths.sumOverLinks(*this);
//   double oldAction = totalAction;
//   spinModelState->flipSpin(ipart);
// #ifdef ENABLE_MPI
//   if (mpi) {
//     spinModelState->broadcastToMPIWorkers(mpi);
//   }
// #endif
//   paths.sumOverLinks(*this);
//   double newAction = totalAction;
//   spinModelState->flipSpin(ipart);
// #ifdef ENABLE_MPI
//   if (mpi) {
//     spinModelState->broadcastToMPIWorkers(mpi);
//   }       
// #endif
// 
//   return newAction-oldAction;
// }

// void SpinChoiceFixedNodeAction::handleLink(const LinkSummable::Vec &start,
//             const LinkSummable::Vec &end, int ipart, int islice, 
// 	    const Paths &paths) {
//   double u=0., utau=0., ulambda=0;
//   FixedNodeAction::Vec fm=0., fp=0.;
//   this->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
//   totalAction += u;
// }


