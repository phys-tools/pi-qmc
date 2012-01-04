#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"

#include "MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class EMARateActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerInterface *sampler;
    Species species1, species2;
    SimulationInfo simInfo;
    static const int npart=2;
    static const int nmoving=2;
    static const int nlevel=6;
    static const int nslice=64;

    void setIdenticalPaths() {
        Beads<NDIM>::Vec electronPosition(0.0);
        Beads<NDIM>::Vec holePosition(1.0);
        Beads<NDIM> sectionBeads = sampler->getSectionBeads();
        Beads<NDIM> movingBeads = sampler->getMovingBeads();
        for (int i = 0; i < nslice; ++i) {
            sectionBeads(0,i) = holePosition;
            sectionBeads(1,i) = electronPosition;
            movingBeads(0,i) = holePosition;
            movingBeads(1,i) = electronPosition;
        }
    }
};

TEST_F(EMARateActionTest, testActionDifferenceForIdenticalPathsIsZero) {
    EMARateAction action(simInfo, species1, species2, 1.0);
    setIdenticalPaths();
    double deltaAction = action.testableGetActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

}

