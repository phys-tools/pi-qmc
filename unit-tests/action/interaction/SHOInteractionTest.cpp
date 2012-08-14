#include <gtest/gtest.h>

#include "action/interaction/SHOInteraction.h"

#include "advancer/MultiLevelSamplerFake.h"
#include "base/Beads.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"

namespace {

class SHOInteractionTest: public testing::Test {
protected:
    typedef Beads<NDIM>::Vec Vec;

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        setupSpecies(species1);
        setupSpecies(species2);
        simInfo.setTau(0.1);
        omega = 1.0;
    }

    virtual void TearDown() {
        delete sampler;
    }

    void setupSpecies(Species *species) {
        species = new Species();
        species->count = 1;
        simInfo.speciesList.push_back(species);
        simInfo.speciesIndex.push_back(species);
    }

    MultiLevelSamplerFake *sampler;
    Species *species1;
    Species *species2;
    SimulationInfo simInfo;
    double omega;
    static const int npart = 2;
    static const int nmoving = 1;
    static const int nlevel = 6;
    static const int nslice = 64;

    void setPaths(Vec newPos1, Vec newPos2, Vec oldPos1, Vec oldPos2) {
        Beads<NDIM> sectionBeads = sampler->getSectionBeads();
        Beads<NDIM> movingBeads = sampler->getMovingBeads();
        for (int i = 0; i < nslice; ++i) {
            sectionBeads(0, i) = oldPos1;
            sectionBeads(1, i) = oldPos2;
            movingBeads(0, i) = newPos1;
            movingBeads(1, i) = newPos2;
        }
    }
};

TEST_F(SHOInteractionTest, getActionDifferenceForIdenticalPathsIsZero) {
    SHOInteraction action(simInfo, omega);
    setPaths(Vec(0.0), Vec(0.0), Vec(0.0), Vec(0.0));
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_DOUBLE_EQ(0.0, deltaAction);
}

}
