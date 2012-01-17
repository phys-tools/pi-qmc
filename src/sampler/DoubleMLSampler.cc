// $Id$
/*  Copyright (C) 2004-2006 John B. Shumway, Jr. and Saad A. Khairallah

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
#include "RandomNumGenerator.h"
#include "DoubleMLSampler.h"
#include "stats/AccRejEstimator.h"
#include "Beads.h"
#include "Action.h"
#include "DoubleAction.h"
#include "Mover.h"
#include "Paths.h"
#include "Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include "RandomPermutationChooser.h"
#include "DoubleSectionChooser.h"
#include "BeadFactory.h"

DoubleMLSampler::DoubleMLSampler(int nmoving, Paths& paths,
  DoubleSectionChooser &sectionChooser,
  ParticleChooser& particleChooser, PermutationChooser& permutationChooser,
  ParticleChooser& particleChooser2, PermutationChooser& permutationChooser2,
  Mover& mover, Action* action, DoubleAction* doubleAction, const bool both,
  const int nrepeat, const BeadFactory &beadFactory, 
  const bool delayedRejection, const double defaultFactor, double newFactor)
  : MultiLevelSampler(nmoving,paths,sectionChooser,particleChooser,
      permutationChooser,mover,action,nrepeat,beadFactory,
      delayedRejection,defaultFactor,newFactor),
    sectionBeads1(&sectionChooser.getBeads(1)),
    sectionBeads2(&sectionChooser.getBeads(2)),
    sectionPermutation1(&sectionChooser.getPermutation(1)),
    sectionPermutation2(&sectionChooser.getPermutation(2)),
    movingBeads1(movingBeads),
    movingBeads2(beadFactory.getNewBeads(nmoving,sectionBeads->getNSlice())), 
    rejectedBeads1(rejectedBeads),
    rejectedBeads2(delayedRejection
      ?beadFactory.getNewBeads(nmoving,sectionBeads->getNSlice()):0),
    doubleAction(doubleAction),
    movingIndex1(movingIndex), movingIndex2(new IArray(nmoving)),
    pMovingIndex2(nmoving),
    doubleSectionChooser(sectionChooser), samplingBoth(both),
    permutation1(nmoving), permutation2(nmoving),
    particleChooser2(particleChooser2),
    permutationChooser2(permutationChooser2),
    nrepeat(nrepeat) {
  for (int i=0; i<nmoving; ++i) (*movingIndex2)(i)=identityIndex(i)=i;

}

DoubleMLSampler::~DoubleMLSampler() {
  delete movingBeads2;  delete rejectedBeads2;
  delete movingIndex2;
  delete &permutationChooser2;
  delete &particleChooser2;
}

void DoubleMLSampler::run() {
  activateSection(1); permutationChooser.init();
  if (samplingBoth) {
    activateSection(2); permutationChooser2.init(); activateSection(1);
  }
  // Select particles to move and the permutation.
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    bool isNewPerm=permutationChooser.choosePermutation();
    permutation1=permutationChooser.getPermutation(); 
    particleChooser.chooseParticles();
    double lnTranProb=permutationChooser.getLnTranProb();
    for (int i=0; i<nmoving; ++i) (*movingIndex1)(i)=particleChooser[i];
    if (samplingBoth) {
      isNewPerm &= permutationChooser2.choosePermutation();
      permutation2=permutationChooser2.getPermutation();
      lnTranProb+=permutationChooser2.getLnTranProb();
      particleChooser2.chooseParticles();
      for (int i=0; i<nmoving; ++i) (*movingIndex2)(i)=particleChooser2[i];
    }
    
    if (isNewPerm) {
  // Copy old coordinate endpoint to the moving coordinate endpoints.
      const int nsectionSlice=movingBeads->getNSlice();
      for (int imoving=0; imoving<nmoving; ++imoving) {
        pMovingIndex(imoving)=(*movingIndex1)(permutation1[imoving]);
        pMovingIndex2(imoving)=(*movingIndex2)(permutation2[imoving]);
      }
      
      sectionBeads1->copySlice(*movingIndex1,0,*movingBeads1,identityIndex,0);
     
   
      sectionBeads1->copySlice(pMovingIndex,nsectionSlice-1,
                               *movingBeads1,identityIndex,nsectionSlice-1);
      if (samplingBoth) {
        sectionBeads2->copySlice(*movingIndex2,0,*movingBeads2,identityIndex,0);
        sectionBeads2->copySlice(pMovingIndex2,nsectionSlice-1,
                                 *movingBeads2,identityIndex,nsectionSlice-1);
      }
      if (tryMove(lnTranProb) && irepeat<nrepeat-1) {
    
        permutationChooser.init();
        if (samplingBoth) {
          activateSection(2); permutationChooser2.init(); activateSection(1); 
        }
      }
    }
  }
}

bool DoubleMLSampler::tryMove(double initialLnTranProb) {
  double oldDeltaAction=initialLnTranProb;
  for (int ilevel=nlevel; ilevel>=0; --ilevel) {
    if (accRejEst) accRejEst->tryingMove(ilevel);
  
    // Make the trial move for this level A2B, and get the acc ratio and transition prob. A2B.
    factor = defaultFactor;
    double lnTranProb=mover.makeMove(*this,ilevel);
    double transA2B = mover.getForwardProb();
    double deltaAction=(action==0)?0:action->getActionDifference(*this,ilevel);
    if (samplingBoth) {
      activateSection(2);
      lnTranProb+=mover.makeMove(*this,ilevel);
      transA2B+=mover.getForwardProb();
      deltaAction+=(action==0)?0:action->getActionDifference(*this,ilevel);
      activateSection(1);
    }
    deltaAction+=doubleAction->getActionDifference(*this,ilevel);
    double piRatioBoA= -deltaAction+oldDeltaAction;    
    double acceptProb=exp(lnTranProb+piRatioBoA);

    // If you want to do DelayedRejection
    if (delayedRejection && ilevel==nlevel-1){
      if (RandomNumGenerator::getRand()>acceptProb) {
       	double accRatioA2B = acceptProb;
	for (int i=0; i<nmoving; ++i) {
	  for (int islice=0; islice<sectionBeads->getNSlice(); ++islice) {
	    (*rejectedBeads1)(i,islice)=(*movingBeads1)(i,islice);
	    (*rejectedBeads2)(i,islice)=(*movingBeads2)(i,islice);
	  }
	}
	
	//Make another trial move for this level A2C, and get the acc ratio.
        factor = newFactor;
        lnTranProb=mover.makeMove(*this,ilevel);
	deltaAction=(action==0)?0:action->getActionDifference(*this,ilevel);
	if (samplingBoth) {
	  activateSection(2);
	  lnTranProb+=mover.makeMove(*this,ilevel);
	  deltaAction+=(action==0)?0:action->getActionDifference(*this,ilevel);
	  activateSection(1);
	}
	deltaAction+=doubleAction->getActionDifference(*this,ilevel);
	double piRatioCoA = -deltaAction+oldDeltaAction;    
	double accRatioA2C = exp(lnTranProb+piRatioCoA);
	
	//nodal crossing A2C
	//	if (deltaAction > 2e90) return false;
	
	//Make the trial move to the rejected state C2B, and get the acc ratio and transition prob. C2B
	factor=defaultFactor;
	lnTranProb=mover.makeDelayedMove(*this,ilevel);
	double transC2B = mover.getForwardProb();
	if (samplingBoth) {
	  activateSection(2);
	  lnTranProb+=mover.makeDelayedMove(*this,ilevel); 
	  transC2B += mover.getForwardProb();
	  activateSection(1);
	}
	double accRatioC2B = exp(lnTranProb - piRatioCoA  + piRatioBoA);
	accRatioC2B=(accRatioC2B>=1)?1:accRatioC2B;
	double accRatioA2B2C = accRatioA2C*exp(transC2B+transA2B)*(1-accRatioC2B)/(1-accRatioA2B);
		deltaAction=deltaAction-transC2B-transA2B-log((1-accRatioC2B)/(1-accRatioA2B));
	if (RandomNumGenerator::getRand()>accRatioA2B2C) return false;
      }
    } else {  
      if (RandomNumGenerator::getRand()>acceptProb) return false;
    }
      oldDeltaAction=deltaAction;
    

    if (accRejEst) accRejEst->moveAccepted(ilevel); 
  
  }

  // Move accepted.
  action->acceptLastMove();

  //if (samplingBoth)  {
  //  activateSection(2); action->acceptLastMove(); activateSection(1);
  // }
  doubleAction->acceptLastMove();
  // Put moved beads in section beads.
  for (int islice=0; islice<sectionBeads->getNSlice(); ++islice) {
    movingBeads1->copySlice(identityIndex,islice,
             *sectionBeads1,*movingIndex1,islice);
    if (samplingBoth)  {
      movingBeads2->copySlice(identityIndex,islice,
               *sectionBeads2,*movingIndex2,islice);
      }
  }
 

  // Append the current permutation to section permutation.
  Permutation temp1(permutation1), temp2(permutation2);
  for (int i=0; i<nmoving; ++i) {
    temp1[i]=(*sectionPermutation1)[(*movingIndex1)(i)];
    if (samplingBoth) temp2[i]=(*sectionPermutation2)[(*movingIndex2)(i)];
  }
  for (int i=0; i<nmoving; ++i) {
    (*sectionPermutation1)[(*movingIndex1)(i)]= temp1[permutation1[i]];
    if (samplingBoth)
      (*sectionPermutation2)[(*movingIndex2)(i)]=temp2[permutation2[i]];
  }


  return true;
}

void DoubleMLSampler::activateSection(const int i) {
  if (i==1) {
    sectionBeads=sectionBeads1;
    sectionPermutation=sectionPermutation1;
    movingBeads=movingBeads1;
    movingIndex=movingIndex1;
    if (delayedRejection) rejectedBeads=rejectedBeads1;
  } else {
    sectionBeads=sectionBeads2;
    sectionPermutation=sectionPermutation2;
    movingBeads=movingBeads2;
    movingIndex=movingIndex2;
    if (delayedRejection) rejectedBeads=rejectedBeads2;
  }
}

AccRejEstimator* 
DoubleMLSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": level " << nlevel << ", moving " << nmoving
           << " " << particleChooser.getName();
  if (samplingBoth) longName << ", both";
  return accRejEst=new AccRejEstimator(longName.str(),nlevel+1);
}
