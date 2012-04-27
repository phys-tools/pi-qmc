#include <gtest/gtest.h>

#include "action/PrimSHOAction.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class PrimSHOActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        species.count = npart;
        simInfo.setTau(0.1);
        a = 0.5;
        b = 0.0;
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species;
    SimulationInfo simInfo;
    double a;
    double b;
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

TEST_F(PrimSHOActionTest, getActionDifferenceForIdenticalPathsIsZero) {
    PrimSHOAction action(a, b, simInfo, NDIM, species);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

TEST_F(PrimSHOActionTest, getActionDifferenceForOneMovedBead) {
    PrimSHOAction action(a, b, simInfo, NDIM, species);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 3.0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    double r2 = dot(position, position);
    double expect = a * r2 * simInfo.getTau();
    ASSERT_FLOAT_EQ(expect, deltaAction);
}

TEST_F(PrimSHOActionTest, getActionDifferenceForOneMovedBeadWithQuadraticTerm) {
    b = 1.0;
    PrimSHOAction action(a, b, simInfo, NDIM, species);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 3.0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    double r2 = dot(position, position);
    double expect = (a * r2 + b * r2 * r2) * simInfo.getTau();
    ASSERT_FLOAT_EQ(expect, deltaAction);
}

}

