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

double SpinChoiceFixedNodeAction::getActionDifference(const Paths &paths,
    int ipart) {
  if (!FixedNodeAction::nodeModel->dependsOnOtherParticles() ) {
    if (ipart==FixedNodeAction::npart-1) return 0;
  }

  int  nsliceOver2=(paths.getNSlice()/2);

  blitz::Range allPart = blitz::Range::all();
  FixedNodeAction::Array prevd1(FixedNodeAction::dist(0,0,allPart));
  FixedNodeAction::Array prevd2(FixedNodeAction::dist(0,1,allPart));
  FixedNodeAction::Array d1(FixedNodeAction::dist(1,0,allPart));
  FixedNodeAction::Array d2(FixedNodeAction::dist(1,1,allPart));

  double deltaAction=0;  
  double oldDet = 0.;
 
  // First check for nodal crossing with new positions 
  // and calculate new nodal action.
  spinModelState->flipSpin(ipart);
#ifdef ENABLE_MPI
  if (mpi) {
    spinModelState->broadcastToMPIWorkers(mpi);
  }
#endif
  int iFirstSlice = paths.getLowestOwnedSlice(true);
  int nSlice = paths.getHighestOwnedSlice(true);
  for (int islice=iFirstSlice; islice<=nSlice; islice++) {
    if (islice<nSlice) {
      for (int i=0; i<FixedNodeAction::npart; ++i) { 
	FixedNodeAction::r1(i)=paths(i,islice);
        FixedNodeAction::r2(i)=paths(i,islice+nsliceOver2);
      }  
    } else { // be sure to grab last slice carefully!
      for (int i=0; i<FixedNodeAction::npart; ++i) { 
        FixedNodeAction::r1(i)=paths(i,islice-1,1);
        FixedNodeAction::r2(i)=paths(i,islice+nsliceOver2-1,1);
      }  
    }
    
    NodeModel::DetWithFlag result = FixedNodeAction::nodeModel->evaluate(FixedNodeAction::r1,FixedNodeAction::r2,0,false);
    if (result.err) return deltaAction = 2e100;

    if (islice==iFirstSlice) oldDet = result.det;

    if (result.det * oldDet < 0.0) return deltaAction=2e100;

    FixedNodeAction::nodeModel->evaluateDistance(FixedNodeAction::r1,FixedNodeAction::r2,0,d1,d2);
    
    if (islice>iFirstSlice) {
      if (FixedNodeAction::useManyBodyDistance) {
        double d02=0., d12=0., d0p2=0., d1p2=0.;
        for (int i=0; i<FixedNodeAction::npart; ++i) {
          d02 += 1./(d1(i)*d1(i));
          d0p2 += 1./(prevd1(i)*prevd1(i));
          d12 += 1./(d2(i)*d2(i));
          d1p2 += 1./(prevd2(i)*prevd2(i));
        }
        deltaAction-=log( (1-exp(-1./sqrt(d02*d0p2)))
                         *(1-exp(-1./sqrt(d12*d1p2))));
      } else {
        for (int i=0; i<FixedNodeAction::npart; ++i) {
          deltaAction -= log((1-exp(-d1(i)*prevd1(i)))
                            *(1-exp(-d2(i)*prevd2(i))));
        }
      }
    }
    for (int i=0; i<FixedNodeAction::npart; ++i) {
      prevd1(i) = d1(i); 
      prevd2(i) = d2(i); 
    }
  }
  // Then calculate old nodal action.
  spinModelState->flipSpin(ipart);
#ifdef ENABLE_MPI
  if (mpi) {
    spinModelState->broadcastToMPIWorkers(mpi);
  }
#endif
  for (int islice=iFirstSlice; islice<=nSlice; islice++) {
    for (int i=0; i<FixedNodeAction::npart; ++i) { 
      FixedNodeAction::r1(i)=paths(i,islice);
      FixedNodeAction::r2(i)=paths(i,islice+nsliceOver2);
    }
    NodeModel::DetWithFlag result= FixedNodeAction::nodeModel->evaluate(FixedNodeAction::r1,FixedNodeAction::r2,0,false);
    if (result.err) return deltaAction=2e100;
    FixedNodeAction::nodeModel->evaluateDistance(FixedNodeAction::r1,FixedNodeAction::r2,0,d1,d2);
    if (islice>iFirstSlice) {
      if (FixedNodeAction::useManyBodyDistance) {
        double d02=0., d12=0., d0p2=0., d1p2=0.;
        for (int i=0; i<FixedNodeAction::npart; ++i) {
          d02 += 1./(d1(i)*d1(i));
          d0p2 += 1./(prevd1(i)*prevd1(i));
          d12 += 1./(d2(i)*d2(i));
          d1p2 += 1./(prevd2(i)*prevd2(i));
        }
        deltaAction+=log( (1-exp(-1./sqrt(d02*d0p2)))
                         *(1-exp(-1./sqrt(d12*d1p2))));
      } else {
        for (int i=0; i<FixedNodeAction::npart; ++i) {
          deltaAction += log((1-exp(-d1(i)*prevd1(i)))
                            *(1-exp(-d2(i)*prevd2(i))));
        }
      }
    }
    for (int i=0; i<FixedNodeAction::npart; ++i) {
      prevd1(i) = d1(i); 
      prevd2(i) = d2(i); 
    }
  }
  
  return deltaAction;
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


