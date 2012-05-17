#include <gtest/gtest.h>

#include "action/GaussianAction.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class GaussianActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        species.count = npart;
        simInfo.setTau(0.1);
        v0 = 50;
        alpha = 0.0;
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species;
    SimulationInfo simInfo;
    double v0;
    double alpha;
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



TEST_F(GaussianActionTest, getActionDifferenceForIdenticalPathsIsZero) {
    GaussianAction action(v0, alpha, simInfo);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

TEST_F(GaussianActionTest, getActionDifferenceForOneMovedBead) {
    alpha = 1.0;
    GaussianAction action(v0, alpha, simInfo);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 1.0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0, deltaAction);
}

}
