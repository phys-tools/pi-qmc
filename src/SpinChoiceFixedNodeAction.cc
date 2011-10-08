//$Id: SpinChoiceFixedNodeAction.cc 383 2011-04-20 16:56:02Z john.shumwayjr $
/*  Copyright (C) 2004-2009, 2011 John B. Shumway, Jr.

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
#include "SpinChoiceFixedNodeAction.h"
#include "SpinModelState.h"
#include "Paths.h"
#include <cstdlib>
#include <blitz/array.h>

SpinChoiceFixedNodeAction::SpinChoiceFixedNodeAction(
  const SimulationInfo &simInfo,
  const Species &species, NodeModel *nodeModel, bool withNodalAction,
  bool useDistDerivative, int maxlevel, bool useManyBodyDistance) 
  : FixedNodeAction(simInfo,species,nodeModel,withNodalAction,
      useDistDerivative,maxlevel,useManyBodyDistance) {
  std::cout << "nSpeciesPart" << nSpeciesPart << std::endl;
  spinModelState = new SpinModelState(nSpeciesPart);
  modelState = spinModelState;
  nodeModel->setSpinModelState(spinModelState);
}

SpinChoiceFixedNodeAction::~SpinChoiceFixedNodeAction() {
  delete spinModelState;
}

void SpinChoiceFixedNodeAction::initCalc(const int nslice, const int firstSlice) {
  actionDifference = 0.;
}

double SpinChoiceFixedNodeAction::getActionDifference(const Paths &paths,
    int ipart) {
//  double oldAction = this->getTotalAction(paths);
//  std::cout<<"Old action = "<<oldAction<<std::endl;
//  std::exit(-1);
//  spinModelState->flipSpin(ipart);
//  double newAction = this->getTotalAction(paths);
//  spinModelState->flipSpin(ipart);
//  return newAction-oldAction;
  paths.sumOverLinks(*this);
  double oldAction = actionDifference;
  spinModelState->flipSpin(ipart);
  paths.sumOverLinks(*this);
  spinModelState->flipSpin(ipart);
  return (actionDifference-oldAction);
}

void SpinChoiceFixedNodeAction::handleLink(const LinkSummable::Vec &start,
            const LinkSummable::Vec &end, int ipart, int islice, 
	    const Paths &paths) {
  double u=0., utau=0., ulambda=0;
  FixedNodeAction::Vec fm=0., fp=0.;
  this->getBeadAction(paths,ipart,islice,u,utau,ulambda,fm,fp);
  actionDifference += u;
}

/*
double SpinChoiceFixedNodeAction::getTotalAction(const Paths& paths) const {
  double totalAction = 0.;
  blitz::Range allPart = blitz::Range::all();
  FixedNodeAction::Array prevd1(dist(0,0,allPart));
  FixedNodeAction::Array prevd2(dist(0,1,allPart));
  FixedNodeAction::Array d1(dist(1,0,allPart));
  FixedNodeAction::Array d2(dist(1,1,allPart));
  int nsliceOver2 = paths.getNSlice()/2;
  int npart = paths.getNPart();
  int nSlice = paths.getNSlice();
  for (int islice=0; islice<nSlice; islice++) {
    for (int i=0; i<npart; ++i) { 
      r1(i)=paths(i,islice);
      r2(i)=paths(i,(islice+nsliceOver2)%nSlice);
    }
    NodeModel::DetWithFlag result= nodeModel->evaluate(r1,r2,0,false);
    if (result.err) return totalAction=2e100;
    nodeModel->evaluateDistance(r1,r2,0,d1,d2);
    if (islice>0) {
      if (FixedNodeAction::useManyBodyDistance) {
        double d02=0., d12=0., d0p2=0., d1p2=0.;
        for (int i=0; i<npart; ++i) {
          d02 += 1./(d1(i)*d1(i));
          d0p2 += 1./(prevd1(i)*prevd1(i));
          d12 += 1./(d2(i)*d2(i));
          d1p2 += 1./(prevd2(i)*prevd2(i));
        }
        totalAction += log( (1-exp(-1./sqrt(d02*d0p2)))
                               *(1-exp(-1./sqrt(d12*d1p2))));
      } else {
        for (int i=0; i<npart; ++i) {
          totalAction += log((1-exp(-d1(i)*prevd1(i)))
                                 *(1-exp(-d2(i)*prevd2(i))));
        }
      }
    }
    for (int i=0; i<npart; ++i) {
      prevd1(i) = d1(i); 
      prevd2(i) = d2(i); 
    }
  }
  std::cout<<totalAction<<std::endl;
  return totalAction;
}
*/
