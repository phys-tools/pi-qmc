#include "EMARateTestBeadPositioner.h"
#include "sampler/test/MultiLevelSamplerFake.h"
#include "Beads.h"

EMARateTestBeadPositioner::EMARateTestBeadPositioner(
        MultiLevelSamplerFake& sampler)
: sampler(sampler) {
}

void EMARateTestBeadPositioner::setIdenticalPaths(double separation) {
    Beads<NDIM>::Vec electronPosition(0.0);
    Beads<NDIM>::Vec holePosition(1.0);
    Beads<NDIM> &sectionBeads(sampler.getSectionBeads());
    Beads<NDIM> &movingBeads(sampler.getMovingBeads());
    for(int islice = 0; islice < sampler.nslice; ++islice){
        sectionBeads(0, islice) = holePosition;
        sectionBeads(1, islice) = electronPosition;
    }
    sampler.copySectionBeadsToMovingBeads();
}

void EMARateTestBeadPositioner::setRecombiningPaths(double separation) {
    Beads<NDIM>::Vec beforePosition(0.0, 0.0, 0.0);
    Beads<NDIM>::Vec afterPosition(1.0, 1.0, 1.0);
    Beads<NDIM> &sectionBeads(sampler.getSectionBeads());
    Beads<NDIM> &movingBeads(sampler.getMovingBeads());
    for (int islice = 0; islice < sampler.nslice; ++islice) {
        sectionBeads(0,islice) = beforePosition;
        sectionBeads(1,islice) = afterPosition;
        if (islice < sampler.nslice/2) {
            movingBeads(0,islice) = beforePosition;
            movingBeads(1,islice) = beforePosition;
        } else if (islice == sampler.nslice/2) {
            movingBeads(0, islice) = beforePosition;
            movingBeads(1, islice) = afterPosition;
        } else {
            movingBeads(0,islice) = afterPosition;
            movingBeads(1,islice) = afterPosition;
        }
    }
}

