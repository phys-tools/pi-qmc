#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include <iostream>

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "EMARateTestBeadPositioner.h"


namespace {

class EMARateActionTest: public testing::Test {
protected:

    virtual void SetUp() {
        simInfo.nslice = 128;
        coefficient = 10.0;
    }

    virtual void TearDown() {
        delete sampler;
        delete positioner;
    }

    MultiLevelSamplerFake *sampler;
    EMARateTestBeadPositioner *positioner;
    Species species1, species2;
    SimulationInfo simInfo;
    static const int npart=2;
    int nmoving;
    static const int nslice=8;
    double coefficient;

    void createFakeSampler(int movingCount) {
        nmoving = movingCount;
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        sampler->firstSliceIndex = simInfo.nslice - nslice / 2;
        positioner = new EMARateTestBeadPositioner(*sampler);
    }

};

TEST_F(EMARateActionTest, testActionDifferenceForIdenticalPathsIsZero) {
    createFakeSampler(2);
    EMARateAction action(simInfo, species1, species2, coefficient);
    positioner->setIdenticalPaths(1.0);
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_DOUBLE_EQ(0.0, deltaAction);
}

TEST_F(EMARateActionTest, testActionDifferenceForRecombiningPaths) {
    createFakeSampler(2);
    EMARateAction action(simInfo, species1, species2, coefficient);
    positioner->setRecombiningPaths(1.0);
    double deltaAction = action.getActionDifference(*sampler, 0);
    double oldAction = -log(1 + coefficient * exp(-3.0));
    double newAction = -log(1 + coefficient * exp(+3.0));
    double expect = newAction - oldAction;
    ASSERT_DOUBLE_EQ(expect, deltaAction);
}

TEST_F(EMARateActionTest, testSamplerWithOnlyElectron) {
    createFakeSampler(1);
    (*sampler->movingIndex)(0) = 1;
    EMARateAction action(simInfo, species1, species2, coefficient);
    positioner->setIdenticalPaths(1.0);
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_DOUBLE_EQ(0.0, deltaAction);
}

}

