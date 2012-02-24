#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "CollectiveSectionSampler.h"
#include "Paths.h"
#include "SectionChooser.h"
#include "Action.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "stats/AccRejEstimator.h"
#include "CollectiveSectionMover.h"


CollectiveSectionSampler::CollectiveSectionSampler(int npart,
        SectionChooser &sectionChooser, Action *action, int nrepeat,
        const BeadFactory &beadFactory, CollectiveSectionMover* mover)
  : nlevel(sectionChooser.getNLevel()), npart(npart), 
    sectionBeads(&sectionChooser.getBeads()), 
    movingBeads(beadFactory.getNewBeads(npart,sectionBeads->getNSlice())),
    movingIndex(new IArray(npart)), identityIndex(npart),
    sectionChooser(sectionChooser), cell(0),
    action(action), nrepeat(nrepeat),
    maximumMovingCount(npart), mover(mover), accRejEst(0) {
  for (int i=0; i<npart; ++i) (*movingIndex)(i) = identityIndex(i) = i;
}

CollectiveSectionSampler::~CollectiveSectionSampler() {
  delete movingBeads;
  delete movingIndex;
}

void CollectiveSectionSampler::run() {

}

bool CollectiveSectionSampler::tryMove() {

}

AccRejEstimator* 
CollectiveSectionSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": moving " << npart;
  return accRejEst=new AccRejEstimator(longName.str(),1);
}
