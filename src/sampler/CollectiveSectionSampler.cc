#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "CollectiveSectionSampler.h"
#include "Paths.h"
#include "SectionChooser.h"
#include "Action.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "SuperCell.h"
#include "stats/AccRejEstimator.h"
#include "CollectiveSectionMover.h"
#include "RandomNumGenerator.h"


CollectiveSectionSampler::CollectiveSectionSampler(int npart,
        SectionChooser &sectionChooser, Action *action, int nrepeat,
        const BeadFactory &beadFactory, CollectiveSectionMover* mover,
	SuperCell *cell)
  : nlevel(sectionChooser.getNLevel()), npart(npart), 
    sectionBeads(&sectionChooser.getBeads()), 
    movingBeads(beadFactory.getNewBeads(npart,sectionBeads->getNSlice())),
    movingIndex(new IArray(npart)), identityIndex(npart),
    sectionChooser(sectionChooser), cell(cell),
    action(action), nrepeat(nrepeat),
    maximumMovingCount(npart), mover(mover), accRejEst(0) {
  for (int i=0; i<npart; ++i) (*movingIndex)(i) = identityIndex(i) = i;
}

CollectiveSectionSampler::~CollectiveSectionSampler() {
  delete movingBeads;
  delete movingIndex;
}

void CollectiveSectionSampler::run() {
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    const int nSectionSlice = movingBeads->getNSlice();
    for (int islice=0; islice<nSectionSlice; ++islice) {
      sectionBeads->copySlice(*movingIndex,islice,
			     *movingBeads,identityIndex,islice);
    }
    tryMove();
  }
}

bool CollectiveSectionSampler::tryMove() {
  if (accRejEst) accRejEst->tryingMove(0);
  double lnTranProb = mover->makeMove(*this,0);
  double deltaAction = (action==0)?0:action->getActionDifference(*this,0);
  double acceptProb = exp(lnTranProb - deltaAction);
  if (RandomNumGenerator::getRand()>acceptProb) return false;
  if (accRejEst) accRejEst->moveAccepted(0);
  action->acceptLastMove();
  int nSectionSlice = sectionBeads->getNSlice();
  for (int islice=0; islice<nSectionSlice; ++islice) {
    movingBeads->copySlice(identityIndex,islice,
			   *sectionBeads,*movingIndex,islice);
  }
  return true;
}

AccRejEstimator* 
CollectiveSectionSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": moving " << npart;
  return accRejEst=new AccRejEstimator(longName.str(),1);
}
