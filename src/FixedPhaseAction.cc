//$Id: FixedPhaseAction.cc,v 1.20 2007/10/03 12:53:56 jshumwa Exp $
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

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
#include "FixedPhaseAction.h"
#include "Beads.h"
#include "DoubleSectionChooser.h"
#include "DoubleMLSampler.h"
#include "PhaseModel.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include <blitz/tinyvec-et.h>

FixedPhaseAction::FixedPhaseAction(const SimulationInfo &simInfo,
  const Species &species, PhaseModel *phaseModel, const int maxlevel) 
  : tau(simInfo.getTau()), npart(simInfo.getNPart()),
    nSpeciesPart(species.count), ifirst(species.ifirst), 
    mass(species.mass),
    r1(npart), r2(npart),
    phi((int)(pow(2,maxlevel)+0.1)+1), newPhi((int)(pow(2,maxlevel)+0.1)+1),
    gradPhi1((int)(pow(2,maxlevel)+0.1)+1), 
    gradPhi2((int)(pow(2,maxlevel)+0.1)+1), 
    newGradPhi1((int)(pow(2,maxlevel)+0.1)+1),
    newGradPhi2((int)(pow(2,maxlevel)+0.1)+1),
    vecPot1((int)(pow(2,maxlevel)+0.1)+1), 
    vecPot2((int)(pow(2,maxlevel)+0.1)+1), 
    newVecPot1((int)(pow(2,maxlevel)+0.1)+1),
    newVecPot2((int)(pow(2,maxlevel)+0.1)+1),
    dmValue((int)(pow(2,maxlevel)+0.1)+1),
    newDMValue((int)(pow(2,maxlevel)+0.1)+1),
    dist((int)(pow(2,maxlevel)+0.1)+1),
    newDist((int)(pow(2,maxlevel+1)+0.1)+1),
    force(npart), phaseModel(phaseModel) {
  for (unsigned int i=0; i<gradPhi1.size(); ++i) {
    gradPhi1[i].resize(npart);
    gradPhi2[i].resize(npart);
    newGradPhi1[i].resize(npart);
    newGradPhi2[i].resize(npart);
    vecPot1[i].resize(npart);
    vecPot2[i].resize(npart);
    newVecPot1[i].resize(npart);
    newVecPot2[i].resize(npart);
  }
}

FixedPhaseAction::~FixedPhaseAction() {
  delete phaseModel;
}

double FixedPhaseAction::getActionDifference(const DoubleMLSampler &sampler,
    int level) {
  // Get ready to move paths.
  double deltaAction=0;
  const Beads<NDIM>& sectionBeads1=sampler.getSectionBeads(1);
  const Beads<NDIM>& sectionBeads2=sampler.getSectionBeads(2);
  const Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  const Beads<NDIM>& movingBeads2=sampler.getMovingBeads(2);
  const IArray& index1=sampler.getMovingIndex(1); 
  const IArray& index2=sampler.getMovingIndex(2); 
  if ( (index1(0)<ifirst||index1(0)>=ifirst+nSpeciesPart))  {
    notMySpecies=true; return 0;
  }
  else notMySpecies=false;
  int nSlice=sectionBeads1.getNSlice();
  const int nMoving=index1.size();
  const int nStride=(int)pow(2,level);
  // Calculate the fixed-phase action.
  for (int islice=nStride; islice<nSlice; islice+=nStride) { 
    // Calculate fixed-phase action for attemped move.
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i,islice);
    for (int i=0; i<nMoving; ++i)
                                r1(index1(i))=movingBeads1(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i,islice);
    if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i)
                                r2(index2(i))=movingBeads2(i,islice);
    phaseModel->evaluate(r1,r2,islice);
    newPhi(islice) = phaseModel->getPhi();
    newGradPhi1[islice] = phaseModel->getGradPhi(1);
    newGradPhi2[islice] = phaseModel->getGradPhi(2);
    newVecPot1[islice] = phaseModel->getVecPot(1);
    newVecPot2[islice] = phaseModel->getVecPot(2);
/*std::cout << "r1: " << r1 << std::endl;
std::cout << "r2: " << r2 << std::endl;
std::cout << "dphi1: " << newGradPhi1[islice] << std::endl;
std::cout << "dphi2: " << newGradPhi2[islice] << std::endl;
std::cout << "a1: " << newVecPot1[islice] << std::endl;
std::cout << "a2: " << newVecPot2[islice] << std::endl;*/
    for (int i=0; i<npart; ++i) r1(i)-=sectionBeads1(i,islice-nStride);
    for (int i=0; i<nMoving; ++i) {
      r1(index1(i))=movingBeads1(i,islice)-movingBeads1(i,islice-nStride);
    }
    for (int i=0; i<npart; ++i) r2(i)-=sectionBeads2(i,islice-nStride);
    if (sampler.isSamplingBoth()) for (int i=0; i<nMoving; ++i) {
      r2(index2(i))=movingBeads2(i,islice)-movingBeads2(i,islice-nStride);
    }
    deltaAction += action(newPhi(islice-nStride),newPhi(islice),r1,r2,
                          newGradPhi1[islice-nStride],newGradPhi1[islice],
                          newGradPhi2[islice-nStride],newGradPhi2[islice],
                          newVecPot1[islice-nStride],newVecPot1[islice],
                          newVecPot2[islice-nStride],newVecPot2[islice],
                          nStride*tau,level==0);
    // Calculate old fixed-phase action.
    for (int i=0; i<npart; ++i) {
      r1(i)=sectionBeads1(i,islice)-sectionBeads1(i,islice-nStride);
      r2(i)=sectionBeads2(i,islice)-sectionBeads2(i,islice-nStride);
    }
    deltaAction -= action(phi(islice-nStride),phi(islice),r1,r2,
                          gradPhi1[islice-nStride],gradPhi1[islice],
                          gradPhi2[islice-nStride],gradPhi2[islice],
                          vecPot1[islice-nStride],vecPot1[islice],
                          vecPot2[islice-nStride],vecPot2[islice],
                          nStride*tau,level==0);
  }
  return deltaAction;
}

double FixedPhaseAction::getTotalAction(const Paths&, const int level) const {
  return 0;
}

double FixedPhaseAction::action(const double phi0, const double phi1,
    const VArray& r1, const VArray& r2, 
    const VArray& gradPhi01, const VArray& gradPhi11, 
    const VArray& gradPhi02, const VArray& gradPhi12, 
    const VArray& a01, const VArray& a11, 
    const VArray& a02, const VArray& a12, 
    const double teff, const bool reject) const {
  // Calculate change in phase; optionally reject if greater than one radian.
  double deltaPhi=phi1-phi0;
  //deltaPhi -= (((int)floor(deltaPhi/PI+1.0))/2)*(2.0*PI);
  deltaPhi -= floor(0.5*(deltaPhi/PI+1.0))*2.0*PI;
  if (reject && fabs(deltaPhi)>1.5) return 1e200;
  // Calculate the distance propagated.
  double d=1e-20;
  for (int i=0; i<npart; ++i) {
    d+=dot(r1(i),r1(i))+dot(r2(i),r2(i));
  }
  d = sqrt(d);
  double dinv=1.0/d;
  // Calculate the partial derivatives parallel to the displacement.
  double parPhi0=0,parPhi1=0,parA0=0,parA1=0;
  for (int i=0; i<npart; ++i) {
    parPhi0 += dot(gradPhi01(i),r1(i))+dot(gradPhi02(i),r2(i));
    parPhi1 += dot(gradPhi11(i),r1(i))+dot(gradPhi12(i),r2(i));
    parA0 += dot(a01(i),r1(i))+dot(a02(i),r2(i));
    parA1 += dot(a11(i),r1(i))+dot(a12(i),r2(i));
  }
  parPhi0*=dinv; parPhi1*=dinv; parA0*=dinv; parA1*=dinv;
  // Calculate the action.
  double action = 1.2*deltaPhi*deltaPhi*dinv*dinv;
  action += (2*parPhi0*parPhi0+2*parPhi1*parPhi1-parPhi0*parPhi1)/15.0;
  action -= 0.2*deltaPhi*(parPhi1+parPhi0)*dinv;
  action -= (parA1+parA0)*deltaPhi*dinv;
  action += (parA1-parA0)*(parPhi1-parPhi0)/6.0;
  action += (parA1*parA1+parA0*parA0+parA1*parA0)/3.0;
  if(action<0) std::cout << "ERROR: FixedPhaseAction already negative" << std::endl;
  for (int i=0; i<npart; ++i) {
    //action += 0.5*dot(gradPhi11(i)-a11(i),gradPhi11(i)-a11(i));
    //action += 0.5*dot(gradPhi01(i)-a01(i),gradPhi01(i)-a01(i));
    //action += 0.5*dot(gradPhi12(i)+a12(i),gradPhi12(i)+a12(i));
    //action += 0.5*dot(gradPhi02(i)+a02(i),gradPhi02(i)+a02(i));
    Vec temp = gradPhi11(i)-a11(i)-gradPhi01(i)+a01(i)
               -r1(i)*dinv*(parPhi1-parA1-parPhi0+parA0);
    action += dot(temp,temp)/3.0;
    Vec temp1 = (gradPhi11(i)-a11(i))-r1(i)*dinv*(parPhi1-parA1);
    Vec temp0 = (gradPhi01(i)-a01(i))-r1(i)*dinv*(parPhi0-parA0);
    action += dot(temp1,temp0);
    temp = gradPhi12(i)-a12(i)-gradPhi02(i)+a02(i)
           -r2(i)*dinv*(parPhi1-parA1-parPhi0+parA0);
    action += dot(temp,temp)/3.0;
    temp1 = (gradPhi12(i)-a12(i))-r2(i)*dinv*(parPhi1-parA1);
    temp0 = (gradPhi02(i)-a02(i))-r2(i)*dinv*(parPhi0-parA0);
    action += dot(temp1,temp0);
  }
  action *= teff/(2*mass);
  if(action<0) std::cout << "ERROR: FixedPhaseAction negative" << std::endl;
  return action;
}

void FixedPhaseAction::getBeadAction(const Paths &paths, int ipart, int islice,
    double &u, double &utau, double &ulambda, Vec &fm, Vec &fp) const {
  int totNSlice=paths.getNSlice();
  fm=0; fp=0; u=utau=ulambda=0;
  // Attribute u and utau to first particle.
  // We only calculate determinants when iPart==0, then store
  // the forces in the force array.
  if (ipart==0) {
    int jslice=(islice+totNSlice/2)%totNSlice;
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice);
    phaseModel->evaluate(r1,r2,1);
    newPhi(1) = phaseModel->getPhi();
    newGradPhi1[1] = phaseModel->getGradPhi(1);
    newGradPhi2[1] = phaseModel->getGradPhi(2);
    newVecPot1[1] = phaseModel->getVecPot(1);
    newVecPot2[1] = phaseModel->getVecPot(2);
    for (int i=0; i<npart; ++i) r1(i)=paths(i,islice,-1);
    for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice,-1);
    phaseModel->evaluate(r1,r2,0);
    newPhi(0) = phaseModel->getPhi();
    newGradPhi1[0] = phaseModel->getGradPhi(1);
    newGradPhi2[0] = phaseModel->getGradPhi(2);
    newVecPot1[0] = phaseModel->getVecPot(1);
    newVecPot2[0] = phaseModel->getVecPot(2);
    utau = 0.5*action(newPhi(0),newPhi(1),r1,r2,
           newGradPhi1[0],newGradPhi1[1], newGradPhi2[0],newGradPhi2[1],
           newVecPot1[0],newVecPot1[1], newVecPot2[0],newVecPot2[1], tau)/tau;
  }
}

void FixedPhaseAction::initialize(const DoubleSectionChooser &chooser) {
  const Beads<NDIM>& sectionBeads1=chooser.getBeads(1);
  const Beads<NDIM>& sectionBeads2=chooser.getBeads(2);
  nslice=sectionBeads1.getNSlice();
  for (int islice=0; islice<nslice; ++islice) {
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads1(i,islice);
    for (int i=0; i<npart; ++i) r2(i)=sectionBeads2(i,islice);
    phaseModel->evaluate(r1,r2,islice);
    phi(islice)=phaseModel->getPhi();
    gradPhi1[islice]= phaseModel->getGradPhi(1);
    gradPhi2[islice]= phaseModel->getGradPhi(2);
    vecPot1[islice]= phaseModel->getVecPot(1);
    vecPot2[islice]= phaseModel->getVecPot(2);
  } 
  newPhi(0)=phi(0);
  newGradPhi1[0]=gradPhi1[0];
  newGradPhi2[0]=gradPhi2[0];
  newVecPot1[0]=vecPot1[0];
  newVecPot2[0]=vecPot2[0];
}

void FixedPhaseAction::acceptLastMove() {
  if (notMySpecies) return;
  for (int i=0; i<nslice; ++i) {
    phi(i)=newPhi(i);
    gradPhi1[i]=newGradPhi1[i];
    gradPhi2[i]=newGradPhi2[i];
    vecPot1[i]=newVecPot1[i];
    vecPot2[i]=newVecPot2[i];
  }
}

const double FixedPhaseAction::PI = 3.14159265358979;
const double FixedPhaseAction::C = 137.0359895;
