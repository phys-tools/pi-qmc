#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "CollectiveSectionSampler.h"
#include "CollectiveSectionMover.h"
#include "SectionChooser.h"
#include "action/Action.h"
#include "base/Paths.h"
#include "base/Beads.h"
#include "base/BeadFactory.h"
#include "stats/AccRejEstimator.h"
#include "util/SuperCell.h"
#include "util/RandomNumGenerator.h"

CollectiveSectionSampler::CollectiveSectionSampler(int npart,
        SectionChooser &sectionChooser, Action *action, int nrepeat,
        const BeadFactory &beadFactory, CollectiveSectionMover* mover,
        SuperCell *cell) :
        nlevel(sectionChooser.getNLevel()),
        npart(npart),
        sectionBeads(&sectionChooser.getBeads()),
        movingBeads(beadFactory.getNewBeads(npart, sectionBeads->getNSlice())),
        movingIndex(new IArray(npart)),
        identityIndex(npart),
        sectionChooser(sectionChooser),
        cell(cell),
        action(action),
        nrepeat(nrepeat),
        maximumMovingCount(npart),
        mover(mover),
        accRejEst(0) {
    for (int i = 0; i < npart; ++i) {
        (*movingIndex)(i) = identityIndex(i) = i;
    }
}

CollectiveSectionSampler::~CollectiveSectionSampler() {
    delete movingBeads;
    delete movingIndex;
}

void CollectiveSectionSampler::initializeMovingBeads() {
    const int nSectionSlice = movingBeads->getNSlice();
    for (int islice = 0; islice < nSectionSlice; ++islice) {
        sectionBeads->copySlice(*movingIndex, islice, *movingBeads,
                identityIndex, islice);
    }
}

void CollectiveSectionSampler::run() {
    for (int irepeat = 0; irepeat < nrepeat; ++irepeat) {
        initializeMovingBeads();
        tryMove();
    }
}

bool CollectiveSectionSampler::tryMove() {
    reportAtempt();
    double lnTranProb = mover->makeMove(*this);
    double deltaAction = action->getActionDifference(*this, 0);
    double acceptProb = exp(lnTranProb - deltaAction);
    bool accepted = RandomNumGenerator::getRand() < acceptProb;
    if (accepted) {
        reportAcceptance();
        action->acceptLastMove();
        updateSectionBeads();
    }
    return accepted;
}

void CollectiveSectionSampler::reportAtempt() const {
    if (accRejEst) {
        accRejEst->tryingMove(0);
    }
}

void CollectiveSectionSampler::updateSectionBeads() {
    int nSectionSlice = sectionBeads->getNSlice();
    for (int islice = 0; islice < nSectionSlice; ++islice) {
        movingBeads->copySlice(identityIndex, islice, *sectionBeads,
                *movingIndex, islice);
    }
}

void CollectiveSectionSampler::reportAcceptance() const {
    if (accRejEst) {
        accRejEst->moveAccepted(0);
    }
}

AccRejEstimator*
CollectiveSectionSampler::getAccRejEstimator(const std::string& name) {
    std::ostringstream longName;
    longName << name << ": moving " << npart;
    return accRejEst = new AccRejEstimator(longName.str(), 1);
}
