// $Id: DoubleMLSampler.cc,v 1.24 2008/11/25 13:11:21 jshumwa Exp $
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
  const int nrepeat, const BeadFactory &beadFactory)
  : MultiLevelSampler(nmoving,paths,sectionChooser,particleChooser,
                      permutationChooser,mover,action,nrepeat,beadFactory),
    sectionBeads1(&sectionChooser.getBeads(1)),
    sectionBeads2(&sectionChooser.getBeads(2)),
    sectionPermutation1(&sectionChooser.getPermutation(1)),
    sectionPermutation2(&sectionChooser.getPermutation(2)),
    movingBeads1(movingBeads),
    movingBeads2(beadFactory.getNewBeads(nmoving,sectionBeads->getNSlice())),
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
  delete movingBeads2;
  delete movingIndex2;
  delete &permutationChooser2;
  delete &particleChooser2;
}

void DoubleMLSampler::run() {
  activateSection(1); permutationChooser.init();
  if (samplingBoth) {
    activateSection(2); permutationChooser2.init(); activateSection(1);
  }
  // Select particles to move and the permuation.
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
    // Do the trial move for this level.
    double lnTranProb=mover.makeMove(*this,ilevel);
    // Evaluate the change in action for this level.
    double deltaAction=(action==0)?0:action->getActionDifference(*this,ilevel);
    // Repeat for other section.
    if (samplingBoth) {
      activateSection(2);
      lnTranProb+=mover.makeMove(*this,ilevel);
      deltaAction+=(action==0)?0:action->getActionDifference(*this,ilevel);
      activateSection(1);
    }
    // Evaluate and add change in DoubleAction. 
    deltaAction+=doubleAction->getActionDifference(*this,ilevel);
    double acceptProb=exp(lnTranProb-deltaAction+oldDeltaAction);
    //std::cout << acceptProb << ", " << lnTranProb << ",  " << ilevel << ", "
    //          << -deltaAction  << ", " << oldDeltaAction << std::endl;
    if (RandomNumGenerator::getRand()>acceptProb) return false;
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
  if (nmoving>0) {
  }
  for (int i=0; i<nmoving; ++i) {
    temp1[i]=(*sectionPermutation1)[(*movingIndex1)(i)];
    if (samplingBoth) temp2[i]=(*sectionPermutation2)[(*movingIndex2)(i)];
  }
  for (int i=0; i<nmoving; ++i) {
    (*sectionPermutation1)[(*movingIndex1)(i)]=temp1[permutation1[i]];
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
  } else {
    sectionBeads=sectionBeads2;
    sectionPermutation=sectionPermutation2;
    movingBeads=movingBeads2;
    movingIndex=movingIndex2;
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
