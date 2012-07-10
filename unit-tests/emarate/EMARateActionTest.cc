#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include <iostream>

#include "advancer/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "EMARateTestBeadPositioner.h"
#include "action/coulomb/CoulombLinkAction.h"


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
    double oldAction = -log(1 + coefficient);
    double kineticDifference = NDIM;
    double newAction = -log(1 + coefficient * exp(kineticDifference));
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

TEST_F(EMARateActionTest, testActionDifferenceWithCoulombAction) {
    createFakeSampler(2);
    EMARateAction action(simInfo, species1, species2, coefficient);
    double epsilon = 1.0;
    int norder = 3;
    action.includeCoulombContribution(epsilon, norder);
    positioner->setRecombiningPaths(1.0);
    double deltaAction = action.getActionDifference(*sampler, 0);

    double q1q2 = -1.0;
    double mu = 0.5;
    double deltaTau = 1.0;
    CoulombLinkAction coulomb(q1q2, epsilon, mu, deltaTau, norder);
    CoulombLinkAction::Vec zero = 0.0;
    CoulombLinkAction::Vec delta = 1.0;
    double oldAction = -log(1 + coefficient);
    double kineticDifference = 0.0;
    double coulombDifference = 2.0 * coulomb.getValue(zero, delta);
    kineticDifference -= -NDIM;
    coulombDifference -= 2.0 * coulomb.getValue(zero, zero);
    double actionDifference = kineticDifference + coulombDifference;
    double newAction = -log(1 + coefficient * exp(actionDifference));
    double expect = newAction - oldAction;

    ASSERT_DOUBLE_EQ(expect, deltaAction);
}

}

