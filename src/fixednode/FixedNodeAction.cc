//$Id$
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
#include "FixedNodeAction.h"
#include "SuperCell.h"
#include <cstdlib>
#include "Beads.h"
#include "sampler/DoubleSectionChooser.h"
#include "sampler/DoubleMLSampler.h"
#include "NodeModel.h"
#include "Paths.h"
#include "SimulationInfo.h"

FixedNodeAction::FixedNodeAction(const SimulationInfo &simInfo,
  const Species &species, NodeModel *nodeModel, bool withNodalAction,
  bool useDistDerivative, int maxlevel, bool useManyBodyDistance) 
  : tau(simInfo.getTau()), npart(simInfo.getNPart()),
    nSpeciesPart(species.count), ifirst(species.ifirst), 
    r1(npart), r2(npart),
    dmValue((int)(pow(2,maxlevel)+0.1)+1),
    newDMValue((int)(pow(2,maxlevel)+0.1)+1),
    dist((int)(pow(2,maxlevel)+0.1)+1,2,npart),
    newDist((int)(pow(2,maxlevel)+0.1)+1,2,npart), force(npart), 
    gradd1(npart,npart), gradd2(npart,npart),
    dim1(npart), dip1(npart), di1(npart), dim2(npart), dip2(npart), di2(npart),
    dotdim1(npart),  dotdi1(npart), dotdim2(npart), dotdi2(npart),
    nodeModel(nodeModel), matrixUpdateObj(nodeModel->getUpdateObj()),
    withNodalAction(withNodalAction),
    useDistDerivative(useDistDerivative),
    nerror(0), useManyBodyDistance(useManyBodyDistance) {
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
      NodeModel::DetWithFlag result= nodeModel->evaluate(r1,r2,islice,false);
      if (result.err) return deltaAction=2e100;
      newDMValue(islice)=result.det;
    }
    if (newDMValue(islice)*dmValue(0)<=1e-200) return deltaAction=2e100;
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
      if (useManyBodyDistance) {
        double d12=0., d22=0., d1p2=0., d2p2=0.;
        double newd12=0., newd22=0., newd1p2=0., newd2p2=0.;
        for (int i=0; i<npart; ++i) {
          d12 += 1./(dist(islice,0,i)*dist(islice,0,i));
          d22 += 1./(dist(islice,1,i)*dist(islice,1,i));
          d1p2 += 1./(dist(islice-1,0,i)*dist(islice-1,0,i));
          d2p2 += 1./(dist(islice-1,1,i)*dist(islice-1,1,i));
          newd12 += 1./(newDist(islice,0,i)*newDist(islice,0,i));
          newd22 += 1./(newDist(islice,1,i)*newDist(islice,1,i));
          newd1p2 += 1./(newDist(islice-1,0,i)*newDist(islice-1,0,i));
          newd2p2 += 1./(newDist(islice-1,1,i)*newDist(islice-1,1,i));
        }
        deltaAction+=log( (1-exp(-1./sqrt(d12*d1p2)))
                         /(1-exp(-1./sqrt(newd12*newd1p2))));
        deltaAction+=log( (1-exp(-1./sqrt(d22*d2p2)))
                         /(1-exp(-1./sqrt(newd22*newd2p2))));
      } else {
        for (int i=0; i<npart; ++i) {
          deltaAction+=log( (1-exp(-dist(islice,0,i)*dist(islice-1,0,i)))
                     /(1-exp(-newDist(islice,0,i)*newDist(islice-1,0,i))) );
          deltaAction+=log( (1-exp(-dist(islice,1,i)*dist(islice-1,1,i)))
                     /(1-exp(-newDist(islice,1,i)*newDist(islice-1,1,i))) );
        }
      }
    }
  }
 
  return deltaAction;
}

double FixedNodeAction::getActionDifference(const Paths &paths,
    const VArray &displacement, int nMoving, const IArray &movingIndex, 
    int iFirstSlice, int nSlice) {
  //Check if species needs fixed node action
  if (!nodeModel->dependsOnOtherParticles() ) {
    for (int i=0; i<nMoving; ++i) {
      if ( movingIndex(i)>=ifirst 
        && movingIndex(i)<ifirst+nSpeciesPart) break;
      if (i==nMoving-1) return 0;
    }
  }

  int  nsliceOver2=(paths.getNSlice()/2);
  const SuperCell& cell=paths.getSuperCell();

  blitz::Range allPart = blitz::Range::all();
  Array prevd1(dist(0,0,allPart));
  Array prevd2(dist(0,1,allPart));
  Array d1(dist(1,0,allPart));
  Array d2(dist(1,1,allPart));

  double deltaAction=0;  
  double oldDet = 0.;
 
  // First check for nodal crossing with new positions 
  // and calculate new nodal action.
  for (int islice=iFirstSlice; islice<=nSlice; islice++) {
//std::cout << "Checking node crossing for slice " << islice << std::endl;
    if (islice<nSlice) {
      for (int i=0; i<npart; ++i) { 
        r1(i)=paths(i,islice);
        r2(i)=paths(i,islice+nsliceOver2);
      }  
    } else { // be sure to grab last slice carefully!
      for (int i=0; i<npart; ++i) { 
        r1(i)=paths(i,islice-1,1);
        r2(i)=paths(i,islice+nsliceOver2-1,1);
      }  
    }
    
    for (int i=0; i<nMoving; ++i) {
      r1(movingIndex(i))+=displacement(i);
      cell.pbc(r1(movingIndex(i)));
      r2(movingIndex(i))+=displacement(i);
      cell.pbc(r2(movingIndex(i)));
    }

    NodeModel::DetWithFlag result = nodeModel->evaluate(r1,r2,0,false);
    if (result.err) return deltaAction = 2e100;

    if (islice==iFirstSlice) oldDet = result.det;

    if (result.det * oldDet < 0.0) return deltaAction=2e100;

    //std::cout << islice << ", " << result.det << std::endl;

    nodeModel->evaluateDistance(r1,r2,0,d1,d2);
    
    if (islice>iFirstSlice) {
      if (useManyBodyDistance) {
        double d02=0., d12=0., d0p2=0., d1p2=0.;
        for (int i=0; i<npart; ++i) {
          d02 += 1./(d1(i)*d1(i));
          d0p2 += 1./(prevd1(i)*prevd1(i));
          d12 += 1./(d2(i)*d2(i));
          d1p2 += 1./(prevd2(i)*prevd2(i));
        }
        deltaAction-=log( (1-exp(-1./sqrt(d02*d0p2)))
                         *(1-exp(-1./sqrt(d12*d1p2))));
      } else {
        for (int i=0; i<npart; ++i) {
          deltaAction -= log((1-exp(-d1(i)*prevd1(i)))
                            *(1-exp(-d2(i)*prevd2(i))));
        }
      }
    }
    for (int i=0; i<npart; ++i) {
      prevd1(i) = d1(i); 
      prevd2(i) = d2(i); 
    }
  }
  // Then calculate old nodal action.
  for (int islice=iFirstSlice; islice<=nSlice; islice++) {
    for (int i=0; i<npart; ++i) { 
      r1(i)=paths(i,islice);
      r2(i)=paths(i,islice+nsliceOver2);
    }
    NodeModel::DetWithFlag result= nodeModel->evaluate(r1,r2,0,false);
    if (result.err) return deltaAction=2e100;
    nodeModel->evaluateDistance(r1,r2,0,d1,d2);
    if (islice>iFirstSlice) {
      if (useManyBodyDistance) {
        double d02=0., d12=0., d0p2=0., d1p2=0.;
        for (int i=0; i<npart; ++i) {
          d02 += 1./(d1(i)*d1(i));
          d0p2 += 1./(prevd1(i)*prevd1(i));
          d12 += 1./(d2(i)*d2(i));
          d1p2 += 1./(prevd2(i)*prevd2(i));
        }
        deltaAction+=log( (1-exp(-1./sqrt(d02*d0p2)))
                         *(1-exp(-1./sqrt(d12*d1p2))));
      } else {
        for (int i=0; i<npart; ++i) {
          deltaAction += log((1-exp(-d1(i)*prevd1(i)))
                            *(1-exp(-d2(i)*prevd2(i))));
        }
      }
    }
    for (int i=0; i<npart; ++i) {
      prevd1(i) = d1(i); 
      prevd2(i) = d2(i); 
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
    NodeModel::DetWithFlag
      detm = nodeModel->evaluate(r1, r2, 0, false);
    nodeModel->evaluateDistance(r1,r2,0,dim1,dim2);
    if (useDistDerivative) {
      nodeModel->evaluateDotDistance(r1,r2,0,dotdim1,dotdim2);
    } else {
      dotdim1=0.; dotdim2=0.;
    }
    // Calculate the action and the gradient of the action.
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice);
    NodeModel::DetWithFlag
      det = nodeModel->evaluate(r1, r2, 0, false);
    nodeModel->evaluateDistance(r1,r2,0,di1,di2);
    if (useDistDerivative) {
      nodeModel->evaluateDotDistance(r1,r2,0,dotdi1,dotdi2);
    } else {
      dotdi1=0.; dotdi2=0.;
    }
    // For now, forces are zero, so we cannot use the virial estimator.
    force=0.0;
    // Calculate u and utau.
    if (useManyBodyDistance) {
      double d12=0., d1m2=0., dotd11=0., dotd1m1=0.;
      for (int i=0; i<npart; ++i) {
        d12 += 1./(di1(i)*di1(i));
        d1m2 += 1./(dim1(i)*dim1(i));
        dotd11 += dotdi1(i)/(di1(i)*di1(i)*di1(i));
        dotd1m1 += dotdim1(i)/(dim1(i)*dim1(i)*dim1(i));
      }
      d12 = 1./sqrt(d12);
      d1m2 = 1./sqrt(d1m2);
      dotd11 *= d12*d12*d12;
      dotd1m1 *= d1m2*d1m2*d1m2;
      double xim1 = d12*d1m2;
      double exim1 = exp(-xim1);
      u -= log( (1-exim1) );
      utau -= (dotd11*d1m2 + d12*dotd1m1)*exim1/(1-exim1); 
    } else {
      for (int jpart=0; jpart<npart; ++jpart) {
        double xim1=dim1(jpart)*di1(jpart);
        double dotxim1=dotdim1(jpart)*di1(jpart)+dim1(jpart)*dotdi1(jpart);
        // Calculate the nodal action.
        u += -log(1-exp(-xim1));
        if (det.err || detm.err || det.det*detm.det<0) u += 1e200;
        utau += -dotxim1*exp(-xim1)/(1-exp(-xim1)); 
      }
    }
  } else {
    // Just check for node crossing.
    int jslice=(islice+totNSlice/2)%totNSlice;
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice,-1);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice,-1);
    NodeModel::DetWithFlag
      detm = nodeModel->evaluate(r1, r2, 0, false);
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice);
    NodeModel::DetWithFlag
      det = nodeModel->evaluate(r1, r2, 0, false);
    if (det.err || detm.err ||det.det*detm.det < 0) u = 1e200; 
  }
  fm=force(ipart); 
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
    NodeModel::DetWithFlag result = nodeModel->evaluate(r1,r2,islice,true);
    dmValue(islice) = result.det;
    if (result.err || dmValue(islice)*dmValue(0)<=1e-200) {
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
