//$Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr.

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
#include "FixedNodeAction.h"

#include <cstdlib>
#include "Beads.h"
#include "DoubleSectionChooser.h"
#include "DoubleMLSampler.h"
#include "NodeModel.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include "DoubleDisplaceMoveSampler.h"

FixedNodeAction::FixedNodeAction(const SimulationInfo &simInfo,
  const Species &species, NodeModel *nodeModel, 
  const bool withNodalAction, const bool useDistDerivative, const int maxlevel) 
  : tau(simInfo.getTau()), npart(simInfo.getNPart()),
    nSpeciesPart(species.count), ifirst(species.ifirst), 
    r1(npart), r2(npart),
    dmValue((int)(pow(2,maxlevel)+0.1)+1),
    newDMValue((int)(pow(2,maxlevel)+0.1)+1),
    dist((int)(pow(2,maxlevel)+0.1)+1,2,npart),
    newDist((int)(pow(2,maxlevel+1)+0.1)+1,2,npart), force(npart), 
    gradd1(npart,npart), gradd2(npart,npart),
    dim1(npart), dip1(npart), di1(npart), dim2(npart), dip2(npart), di2(npart),
    dotdim1(npart),  dotdi1(npart), dotdim2(npart), dotdi2(npart),
    nodeModel(nodeModel), matrixUpdateObj(nodeModel->getUpdateObj()),
    withNodalAction(withNodalAction),
    useDistDerivative(useDistDerivative),
    nerror(0) {
  std::cout << "FixedNodeAction" << std::endl;
}

FixedNodeAction::~FixedNodeAction() {
  delete nodeModel;
}

double FixedNodeAction::getActionDifference(const DoubleMLSampler &sampler,
    int level) {
  // Get ready to move paths.
  double deltaAction=0;
  const Beads<NDIM>& sectionBeads1=sampler.getSectionBeads(1);
  const Beads<NDIM>& sectionBeads2=sampler.getSectionBeads(2);
  const Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  const Beads<NDIM>& movingBeads2=sampler.getMovingBeads(2);
  const IArray& index1=sampler.getMovingIndex(1); 
  const IArray& index2=sampler.getMovingIndex(2); 
  const int nMoving=index1.size();
  if (!nodeModel->dependsOnOtherParticles() ) {
    for (int i=0; i<nMoving; ++i) {
      if ( (index1(i)>=ifirst && index1(i)<ifirst+nSpeciesPart)) break;
      if (i==nMoving-1) {notMySpecies=true; return 0;}
    }
  }
  notMySpecies=false;
  const int nSlice=sectionBeads1.getNSlice();
  const int nStride=(int)pow(2,level+1);
  // First check for node crossing.
  for (int islice=nStride/2; islice<nSlice; islice+=nStride) {
    if (matrixUpdateObj) {
      newDMValue(islice)=matrixUpdateObj->evaluateChange(sampler, islice)
                        *dmValue(islice);
    } else {
      for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i,islice);
      for (int i=0; i<nMoving; ++i) r1(index1(i))=movingBeads1(i,islice);
      for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i,islice);
      if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                    r2(index2(i))=movingBeads2(i,islice);
      newDMValue(islice)=nodeModel->evaluate(r1,r2,islice);
    }
    if (newDMValue(islice)*dmValue(0)<=0) return deltaAction=2e100;
  } 
  // Calculate the nodal action if level=0;
  if (withNodalAction && level==0) {
    blitz::Range allPart = blitz::Range::all();
    for (int islice=1;islice<nSlice; ++islice) { 
      if (matrixUpdateObj) {
        matrixUpdateObj->evaluateNewInverse(islice);
        Array d1(newDist(islice,0,allPart));
        Array d2(newDist(islice,1,allPart));
        for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i,islice);
        for (int i=0; i<nMoving; ++i)
                                    r1(index1(i))=movingBeads1(i,islice);
        for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i,islice);
        if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                    r2(index2(i))=movingBeads2(i,islice);
        matrixUpdateObj->evaluateNewDistance(r1,r2,islice,d1,d2);
      } else {
        for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i,islice);
        for (int i=0; i<nMoving; ++i)
                                    r1(index1(i))=movingBeads1(i,islice);
        for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i,islice);
        if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                    r2(index2(i))=movingBeads2(i,islice);
        Array d1(newDist(islice,0,allPart));
        Array d2(newDist(islice,1,allPart));
        nodeModel->evaluateDistance(r1,r2,islice,d1,d2);
      }
      for (int i=0; i<npart; ++i) {
//std::cout << dist(islice,0,i) << ", " << dist(islice-1,0,i)  << std::endl;
//std::cout << newDist(islice,0,i) << ", " << newDist(islice-1,0,i)  << std::endl;
//std::cout << dist(islice,1,i) << ", " << dist(islice-1,1,i)  << std::endl;
//std::cout << newDist(islice,1,i) << ", " << newDist(islice-1,1,i)  << std::endl;
        deltaAction+=log( (1-exp(-dist(islice,0,i)*dist(islice-1,0,i)))
                   /(1-exp(-newDist(islice,0,i)*newDist(islice-1,0,i))) );
        deltaAction+=log( (1-exp(-dist(islice,1,i)*dist(islice-1,1,i)))
                   /(1-exp(-newDist(islice,1,i)*newDist(islice-1,1,i))) );
      }
    }
  }
 
  return deltaAction;
}


///////////// displace
double FixedNodeAction::getActionDifference(const DoubleDisplaceMoveSampler &sampler,
   const int nMoving) {
  // Get ready to move paths.
  double deltaAction=0;
  const Beads<NDIM>& pathsBeads1=sampler.getPathsBeads(1);
  const Beads<NDIM>& pathsBeads2=sampler.getPathsBeads(2);
  const Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  const Beads<NDIM>& movingBeads2=sampler.getMovingBeads(2);
  const IArray& index=sampler.getMovingIndex(); 

  if (!nodeModel->dependsOnOtherParticles() ) {
    for (int i=0; i<nMoving; ++i) {
      if ( (index(i)>=ifirst && index(i)<ifirst+nSpeciesPart)) break;
      if (i==nMoving-1) {notMySpecies=true; return 0;}
    }
  }
  notMySpecies=false;
  const int nSlice=pathsBeads1.getNSlice();
  const int nStride=2;// level is 0;(int)pow(2,level+1);
  // First check for node crossing.
  for (int islice=nStride/2; islice<nSlice; islice+=nStride) {
    /*  if (matrixUpdateObj) {
      newDMValue(islice)=matrixUpdateObj->evaluateChange(sampler, islice)
                        *dmValue(islice);
    } else {
    */
      for (int i=0; i<npart; ++i) r1(i)=pathsBeads1(i,islice);
      for (int i=0; i<nMoving; ++i) r1(index(i))=movingBeads1(i,islice);
      for (int i=0; i<npart; ++i) r2(i)=pathsBeads2(i,islice);
      for (int i=0; i<nMoving; ++i) r2(index(i))=movingBeads2(i,islice);
      newDMValue(islice)=nodeModel->evaluate(r1,r2,islice);
      //} matrixupdate
    if (newDMValue(islice)*dmValue(0)<=0) return deltaAction=2e100;
  } 
  // Calculate the nodal action if level=0;
  if (withNodalAction) {
    blitz::Range allPart = blitz::Range::all();
    for (int islice=1;islice<nSlice; ++islice) { 
      /*if (matrixUpdateObj) { not implemented yet
        matrixUpdateObj->evaluateNewInverse(islice);
        Array d1(newDist(islice,0,allPart));
        Array d2(newDist(islice,1,allPart));
        for (int i=0; i<npart; ++i) r1(i)=pathsBeads1(i,islice);
        for (int i=0; i<nMoving; ++i)
                                    r1(index1(i))=movingBeads1(i,islice);
        for (int i=0; i<npart; ++i) r2(i)=pathsBeads2(i,islice);
        if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                    r2(index2(i))=movingBeads2(i,islice);
        matrixUpdateObj->evaluateNewDistance(r1,r2,islice,d1,d2);
      } else {
      */
        for (int i=0; i<npart; ++i) r1(i)=pathsBeads1(i,islice);
        for (int i=0; i<nMoving; ++i) r1(index(i))=movingBeads1(i,islice);
        for (int i=0; i<npart; ++i) r2(i)=pathsBeads2(i,islice);
	for (int i=0; i<nMoving; ++i) r2(index(i))=movingBeads2(i,islice);
        Array d1(newDist(islice,0,allPart));
        Array d2(newDist(islice,1,allPart));
        nodeModel->evaluateDistance(r1,r2,islice,d1,d2);
	//} matrixUpdateObj
      for (int i=0; i<npart; ++i) {
//std::cout << dist(islice,0,i) << ", " << dist(islice-1,0,i)  << std::endl;
//std::cout << newDist(islice,0,i) << ", " << newDist(islice-1,0,i)  << std::endl;
//std::cout << dist(islice,1,i) << ", " << dist(islice-1,1,i)  << std::endl;
//std::cout << newDist(islice,1,i) << ", " << newDist(islice-1,1,i)  << std::endl;
        deltaAction+=log( (1-exp(-dist(islice,0,i)*dist(islice-1,0,i)))
                   /(1-exp(-newDist(islice,0,i)*newDist(islice-1,0,i))) );
        deltaAction+=log( (1-exp(-dist(islice,1,i)*dist(islice-1,1,i)))
                   /(1-exp(-newDist(islice,1,i)*newDist(islice-1,1,i))) );
      }
    }
  }   
 


  return deltaAction;
}





double FixedNodeAction::getTotalAction(const Paths&, const int level) const {
  return 0;
}

void FixedNodeAction::getBeadAction(const Paths &paths, int ipart, int islice,
    double &u, double &utau, double &ulambda, Vec &fm, Vec &fp) const {
  int totNSlice=paths.getNSlice();
  fm=0; fp=0; u=utau=ulambda=0;
  // Attribute u and utau to first particle.
  // We only calculate determinants when iPart==0, then store
  // the forces in the force array.
  if (withNodalAction && ipart==0) {
    // Calculate the action and the gradient of the action.
    // Calculate d_i-1
    int jslice=(islice+totNSlice/2)%totNSlice;
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice,-1);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice,-1);
    nodeModel->evaluate(r1, r2, 0);
    nodeModel->evaluateDistance(r1,r2,0,dim1,dim2);
    if (useDistDerivative) {
      nodeModel->evaluateDotDistance(r1,r2,0,dotdim1,dotdim2);
    } else {
      dotdim1=0.; dotdim2=0.;
    }
    // Calculate d_i+1 (NOTE: commented out calculation of force)
//  for (int i=0; i<npart; ++i) r1(i)=paths(i,islice,+1);
//  for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice,+1);
//  nodeModel->evaluate(r1, r2, 0);
//  nodeModel->evaluateDistance(r1,r2,0,dip1,dip2);
    // Calculate the action and the gradient of the action.
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice);
    nodeModel->evaluate(r1, r2, 0);
    nodeModel->evaluateDistance(r1,r2,0,di1,di2);
    if (useDistDerivative) {
      nodeModel->evaluateDotDistance(r1,r2,0,dotdi1,dotdi2);
    } else {
      dotdi1=0.; dotdi2=0.;
    }
    // Now calculate gradient of log of distance to node.
//  nodeModel->evaluateGradLogDist(r1,r2,0,gradd1,gradd2,di1,di2);
    // And calculate the time derivative.
//  nodeModel->evaluateDotDistance(r1,r2,0,dotdi1,dotdi2);
    // Now calulate forces.
    force=0.0;
    for (int jpart=0; jpart<npart; ++jpart) {
//    double xip1=dip1(jpart)*di1(jpart);
      double xim1=dim1(jpart)*di1(jpart);
      double dotxim1=dotdim1(jpart)*di1(jpart)+dim1(jpart)*dotdi1(jpart);
//    double xip2=dip2(jpart)*di2(jpart);
//    double xim2=dim2(jpart)*di2(jpart);
//    for (int ipart=0; ipart<npart; ++ipart) {
//      force(ipart)+=gradd1(ipart,jpart)*(xip1*exp(-xip1)/(1-exp(-xip1))
//                                        +xim1*exp(-xim1)/(1-exp(-xim1)));
//      force(ipart)+=gradd2(ipart,jpart)*(xip2*exp(-xip2)/(1-exp(-xip2))
//                                        +xim2*exp(-xim2)/(1-exp(-xim2)));
//    }
      // Calculate the nodal action.
      u += -log(1-exp(-xim1));
      utau += xim1*exp(-xim1)/(tau*(1-exp(-xim1)));
      utau += -dotxim1*exp(-xim1)/(1-exp(-xim1)); 
//std :: cout << "FNA :: "<<jpart<<". xim1 "<<xim1<<". dotxim1   "<<dotxim1<<". utau  "<<utau<<". u  "<<u<<". tau  "<<tau<<std ::endl;
    }
  }
  fm=force(ipart); 
  //std :: cout << "FNA :: "<<ipart<<" "<<islice<<"  "<<utau<<"  "<<u<<std ::endl;
}

void FixedNodeAction::initialize(const DoubleSectionChooser &chooser) {
  const Beads<NDIM>& sectionBeads1=chooser.getBeads(1);
  const Beads<NDIM>& sectionBeads2=chooser.getBeads(2);
  nslice=sectionBeads1.getNSlice();
  blitz::Range allPart = blitz::Range::all();
  blitz::Range both = blitz::Range::all();
  for (int islice=0; islice<nslice; ++islice) {  
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i,islice);
    dmValue(islice)=nodeModel->evaluate(r1,r2,islice);
    if (dmValue(islice)*dmValue(0)<=0.0) {
      std::cout << "ERROR - crossed node" << islice << std::endl;
      nerror++;
      if (nerror>1000) {
        std::cout << "too many node crossings, exiting" << std::endl;
        std::exit(-1);
      }
    }
    if (withNodalAction) {
      Array d1(dist(islice,0,allPart)), d2(dist(islice,1,allPart));
      nodeModel->evaluateDistance(r1,r2,islice,d1,d2);
    }
  } 
  newDMValue(0)=dmValue(0); newDist(0,both,allPart)=dist(0,both,allPart);
}


// displacemove
void FixedNodeAction::initialize(const DoubleDisplaceMoveSampler &sampler) {
  const Beads<NDIM>& pathsBeads1=sampler.getPathsBeads(1);
  const Beads<NDIM>& pathsBeads2=sampler.getPathsBeads(2);
  nslice=pathsBeads1.getNSlice();
  blitz::Range allPart = blitz::Range::all(); 
  blitz::Range both = blitz::Range::all();
  for (int islice=0; islice<nslice; ++islice) {
    for (int i=0; i<npart; ++i) r1(i)=pathsBeads1(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=pathsBeads2(i,islice);
    dmValue(islice)=nodeModel->evaluate(r1,r2,islice);
    if (dmValue(islice)*dmValue(0)<=0.0) {
      std::cout << "ERROR - crossed node" << islice << std::endl;
      nerror++;
      if (nerror>1000) {
        std::cout << "too many node crossings, exiting" << std::endl;
        std::exit(-1);
      }
    }
    if (withNodalAction) {
      Array d1(dist(islice,0,allPart)), d2(dist(islice,1,allPart));
      nodeModel->evaluateDistance(r1,r2,islice,d1,d2);
    }
  } 
  newDMValue(0)=dmValue(0); newDist(0,both,allPart)=dist(0,both,allPart);
}

void FixedNodeAction::acceptLastMove() {
  if (notMySpecies) return;
  blitz::Range allPart = blitz::Range::all();
  blitz::Range both = blitz::Range::all();
  for (int i=0; i<nslice; ++i) {
    dmValue(i)=newDMValue(i); dist(i,both,allPart)=newDist(i,both,allPart);
  }
  if (matrixUpdateObj) matrixUpdateObj->acceptLastMove(nslice);
}
