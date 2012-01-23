#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "DoubleCollectiveSectionSampler.h"



DoubleCollectiveSectionSampler::~DoubleCollectiveSectionSampler() {
}

DoubleCollectiveSectionSampler::DoubleCollectiveSectionSampler(
        int maximumMovingCount,
        Paths *paths,
        DoubleSectionChooser *sectionChooser, Action *action,
        DoubleAction *doubleAction, const BeadFactory *beadFactory)
    :   maximumMovingCount(maximumMovingCount) {
}

void DoubleCollectiveSectionSampler::run() {

}

const Beads<NDIM> & DoubleCollectiveSectionSampler::getSectionBeads() const {
}

const Beads<NDIM> & DoubleCollectiveSectionSampler::getMovingBeads() const {
}

Beads<NDIM> & DoubleCollectiveSectionSampler::getSectionBeads() {
}

Beads<NDIM> & DoubleCollectiveSectionSampler::getMovingBeads() {
}

const SectionSamplerInterface::IArray&
DoubleCollectiveSectionSampler::getMovingIndex() const {
}

int DoubleCollectiveSectionSampler::getFirstSliceIndex() const {
}

const SuperCell & DoubleCollectiveSectionSampler::getSuperCell() const {
}
