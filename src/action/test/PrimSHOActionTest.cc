#include <gtest/gtest.h>

#include "action/PrimSHOAction.h"

//#include "MultiLevelSamplerFake.h"
#include "sampler/SectionSamplerInterface.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class PrimSHOActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = 0;
//        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
    }

    virtual void TearDown() {
        delete sampler;
    }

    SectionSamplerInterface *sampler;
    Species species;
    SimulationInfo simInfo;
    static const int npart=1;
    static const int nmoving=1;
    static const int nlevel=6;
    static const int nslice=64;

    void setIdenticalPaths() {
        Beads<NDIM>::Vec position(0.0);
        Beads<NDIM> sectionBeads = sampler->getSectionBeads();
        Beads<NDIM> movingBeads = sampler->getMovingBeads();
        for (int i = 0; i < nslice; ++i) {
            sectionBeads(0,i) = position;
            movingBeads(0,i) = position;
        }
    }
};

TEST_F(PrimSHOActionTest, testActionDifferenceForIdenticalPathsIsZero) {
    PrimSHOAction action(1.0, 0.0, simInfo, NDIM, species);
    setIdenticalPaths();
//    double deltaAction = action.testableGetActionDifference(*sampler, 0);
//    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

}

