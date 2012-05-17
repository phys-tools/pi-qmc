#include <gtest/gtest.h>

#include "action/GateAction.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"
#include <math.h>


namespace {

class GateActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        species.count = npart;
        simInfo.setTau(0.1);
        Gvolt = 1;
        sx = 1;
        sy = 1;
        xwidth = 2;
        ywidth = 2;
        xoffset = 1;
        yoffset = 1;
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species;
    SimulationInfo simInfo;
    double Gvolt;
    double sx;
    double sy;
    double xwidth;
    double ywidth;
    double xoffset;
    double yoffset;
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

TEST_F(GateActionTest, getActionDifferenceForIdenticalPathsIsZero) {
    GateAction action(simInfo, Gvolt, sx, sy, xwidth, ywidth, xoffset, yoffset, species);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

TEST_F(GateActionTest, getActionDifferenceForOneMovedBead) {
    GateAction action(simInfo, Gvolt, sx, sy, xwidth, ywidth, xoffset, yoffset, species);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    double normalconst = (tanh(1)-tanh(-1))*(tanh(1)-tanh(-1));
    double expect = ((tanh(1)-tanh(-1))*tanh(2)-tanh(2)*tanh(2))*0.1/normalconst;
    ASSERT_FLOAT_EQ(expect, deltaAction);
}

}
