#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "DoubleCollectiveSectionSampler.h"
#include "Paths.h"
#include "DoubleSectionChooser.h"
#include "Action.h"
#include "DoubleAction.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "stats/AccRejEstimator.h"
#include "CollectiveSectionMover.h"
#include "RandomNumGenerator.h"

/*
   TO DO: Add samplingBoth
*/

DoubleCollectiveSectionSampler::DoubleCollectiveSectionSampler(int npart,
        DoubleSectionChooser& sectionChooser, Action *action,
        DoubleAction *doubleAction, int nrepeat, const BeadFactory &beadFactory,
	CollectiveSectionMover *mover, bool both)
    : CollectiveSectionSampler(npart, sectionChooser,
	            action, nrepeat, beadFactory, mover),
      sectionBeads1(&sectionChooser.getBeads(1)),
      sectionBeads2(&sectionChooser.getBeads(2)),
      movingBeads1(movingBeads),
      movingBeads2(beadFactory.getNewBeads(npart,sectionBeads->getNSlice())),
      movingIndex1(movingIndex), movingIndex2(new IArray(npart)),
      doubleAction(doubleAction), samplingBoth(both) {
  for (int i=0; i<npart; ++i) (*movingIndex2)(i) = i;
}

DoubleCollectiveSectionSampler::~DoubleCollectiveSectionSampler() {
//  delete sectionBeads2;
  delete movingBeads2;
  delete movingIndex2;
}

void DoubleCollectiveSectionSampler::run() {
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    const int nsectionSlice = movingBeads->getNSlice();
    sectionBeads1->copySlice(*movingIndex1, 0,
	           *movingBeads1, identityIndex, 0);
    sectionBeads1->copySlice(*movingIndex1, nsectionSlice-1,
	           *movingBeads1, identityIndex, nsectionSlice-1);
    if (samplingBoth) {
      sectionBeads2->copySlice(*movingIndex2, 0,
	             *movingBeads2, identityIndex, 0);
      sectionBeads2->copySlice(*movingIndex2, nsectionSlice-1,
	             *movingBeads2, identityIndex, nsectionSlice-1);
    }
    tryMove();
  }
}

bool DoubleCollectiveSectionSampler::tryMove() {
  if (accRejEst) accRejEst->tryingMove(0);
  double lnTranProb = mover->makeMove(*this,0);
  double deltaAction = (action==0)?0:action->getActionDifference(*this,0);
  if (samplingBoth) {
    sectionBeads = sectionBeads2;
    movingBeads = movingBeads2;
    movingIndex = movingIndex2;
    lnTranProb += mover->makeMove(*this,0);
    deltaAction += (action==0)?0:action->getActionDifference(*this,0);
    sectionBeads = sectionBeads1;
    movingBeads = movingBeads1;
    movingIndex = movingIndex1;
  }
  deltaAction += doubleAction->getActionDifference(*this,0);
  double acceptProb = exp(lnTranProb - deltaAction);
  if (RandomNumGenerator::getRand()>acceptProb) return false;
  if (accRejEst) accRejEst->moveAccepted(0);
  action->acceptLastMove();
  doubleAction->acceptLastMove();
  int nSectionSlice = sectionBeads->getNSlice();
  for (int islice=0; islice<nSectionSlice; ++islice) {
    movingBeads1->copySlice(identityIndex,islice,
             *sectionBeads1,*movingIndex1,islice);
    if (samplingBoth)  {
      movingBeads2->copySlice(identityIndex,islice,
               *sectionBeads2,*movingIndex2,islice);
    }
  }
//  movingBeads1->copySlice(identityIndex,nSectionSlice,
//                          *sectionBeads1,*movingIndex1,nSectionSlice);
//  if (samplingBoth) {
//    movingBeads2->copySlice(identityIndex,nSectionSlice,
//	                    *sectionBeads2,*movingIndex2,nSectionSlice);
//  }
  return true;
}

AccRejEstimator* 
DoubleCollectiveSectionSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": moving " << npart;
//  if (samplingBoth) longName << ", both";
  return accRejEst=new AccRejEstimator(longName.str(),1);
}

