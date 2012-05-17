#include <gtest/gtest.h>
#include <action/Action.h>
#include "action/GaussianDotAction.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class GaussianDotActionTest: public ::testing::Test {
public:
    typedef blitz::TinyVector<double,NDIM> Vec;
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        species.count = npart;
        simInfo.setTau(0.1);
        v0 = 1;
        alpha = 0.0;
        center = (0, 0, 0);
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species;
    SimulationInfo simInfo;
    double v0;
    double alpha;
    Vec center;
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



TEST_F(GaussianDotActionTest, getActionDifferenceForIdenticalPathsIsZero) {
    GaussianDotAction action(v0, alpha, center, simInfo);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

TEST_F(GaussianDotActionTest, getActionDifferenceForOneMovedBead) {
    alpha = 1.0;
    GaussianDotAction action(v0, alpha, center, simInfo);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 3.0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0, deltaAction);
}

}
