#include "EMARateTestBeadPositioner.h"
#include "sampler/test/MultiLevelSamplerFake.h"
#include "Beads.h"

EMARateTestBeadPositioner::EMARateTestBeadPositioner(
        MultiLevelSamplerFake& sampler)
: sampler(sampler) {
}

void EMARateTestBeadPositioner::setIdenticalPaths(double separation) {
    Beads<NDIM>::Vec electronPosition(0.0);
    Beads<NDIM>::Vec holePosition(separation, separation, separation);
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
    Beads<NDIM>::Vec afterPosition(separation, separation, separation);
    Beads<NDIM> &sectionBeads(sampler.getSectionBeads());
    Beads<NDIM> &movingBeads(sampler.getMovingBeads());
    for (int islice = 0; islice < sampler.nslice; ++islice) {
        double x = double(islice) / (sampler.nslice - 1);
        Beads<NDIM>::Vec interpolation =
                (1.0 - x) * beforePosition + x * afterPosition;
        sectionBeads(0,islice) = interpolation;
        sectionBeads(1,islice) = interpolation;
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

