#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include "SimulationInfo.h"
#include "emarate/EMARateMover.h"
#include "DeterministicEMARateMover.h"
#include "sampler/test/MultiLevelSamplerFake.h"
#include "EMARateTestBeadPositioner.h"
#include "Species.h"

namespace {

class EMARateMoverTest: public testing::Test {
protected:


    virtual void SetUp() {
        nmoving = 2;
        separation = 1.0;
        setupSimulationInfo();
        coefficient = 10.0;
        mover = new DeterministicEMARateMover(simInfo.tau, 1.0, 1.0, maxlevel, coefficient);
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        positioner = new EMARateTestBeadPositioner(*sampler);
    }

    virtual void TearDown() {
        delete mover;
        delete sampler;
        delete positioner;
    }

    void setupSimulationInfo() {
        species1 = new Species("h", 1, 1.0, 0.0, 0, 0);
        species2 = new Species("e", 1, 1.0, 0.0, 0, 0);
        simInfo.tau = 0.1;
        simInfo.nslice = 128;
        simInfo.npart = 2;
        simInfo.speciesList.resize(simInfo.npart);
        simInfo.speciesList[0] = species1;
        simInfo.speciesList[1] = species2;
        simInfo.speciesIndex.resize(simInfo.npart);
        simInfo.speciesIndex[0] = species1;
        simInfo.speciesIndex[1] = species2;
    }

    DeterministicEMARateMover *mover;
    MultiLevelSamplerFake *sampler;
    EMARateTestBeadPositioner *positioner;

    Species *species1;
    Species *species2;
    double coefficient;
    static const int npart = 2;
    int nmoving;
    static const int maxlevel = 5;
    static const int nslice = 1 << maxlevel;
    SimulationInfo simInfo;
    double separation;

    double calculateRadiatingProbability(double delta2) {
        double reducedMass = 0.5;
        double tau = simInfo.tau * nslice / 2;
        double sigma2 = tau / reducedMass;
        return exp(-0.5 * delta2 / sigma2);
    }

    double calculateDiagonalProbability(double delta2) {
        double mass = 1.0;
        double tau = simInfo.tau * nslice;
        double sigma2 = tau / mass;
        return exp(-0.5 * delta2 / sigma2);
    }
};

TEST_F(EMARateMoverTest, testParticleChooser) {
    mover->chooseParticles();
    int holeIndex = (*mover)[0];
    int electronIndex = (*mover)[1];
    ASSERT_EQ(0, holeIndex);
    ASSERT_EQ(1, electronIndex);
}

TEST_F(EMARateMoverTest, testPermutationChooser) {
    mover->choosePermutation();
    ASSERT_TRUE(mover->getPermutation().isIdentity());
}

TEST_F(EMARateMoverTest, testRadiatingProbabilityForIdenticalPaths) {
    positioner->setIdenticalPaths(separation);
    double probability =
            mover->calculateRadiatingProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell(), nslice);
    double delta2 = 3 * separation * separation;
    double expect = pow(calculateRadiatingProbability(delta2), 2);
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testRadiatingProbabilityForRecombiningPaths) {
    positioner->setRecombiningPaths(separation);
    double probability =
            mover->calculateRadiatingProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell(), nslice);
    double expect = 1.0;
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testDiagonalProbabilityForIdenticalPaths) {
    positioner->setIdenticalPaths(separation);
    double probability =
            mover->calculateDiagonalProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell(), nslice);
    double expect = 1.0;
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testDiagonalProbabilityForRecombiningPaths) {
    positioner->setRecombiningPaths(separation);
    double probability =
            mover->calculateDiagonalProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell(), nslice);
    double delta2 = 3 * separation * separation;
    double expect = pow(calculateDiagonalProbability(delta2), 2);
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testDiagonalDecisionForRecombiningPaths) {
    positioner->setRecombiningPaths(separation);
    double delta2 = 3 * separation * separation;
    double diagonalProbability =
            pow(calculateDiagonalProbability(delta2), 2);
    double radiatingProbability = 1.0;
    double threshhold = diagonalProbability
            / (diagonalProbability + coefficient * radiatingProbability);
    mover->nextRandomNumber = threshhold - 1e-6;
    bool isRadiating = mover->chooseDiagonalOrRadiating(
            sampler->getMovingBeads(),
            nslice, sampler->getSuperCell(), nslice);
    ASSERT_FALSE(isRadiating);
}

TEST_F(EMARateMoverTest, testRadiatinglDecisionForRecombiningPaths) {
    positioner->setRecombiningPaths(separation);
    double delta2 = 3 * separation * separation;
    double diagonalProbability =
            pow(calculateDiagonalProbability(delta2), 2);
    double radiatingProbability = 1.0;
    double threshhold = diagonalProbability
            / (diagonalProbability + coefficient * radiatingProbability);
    mover->nextRandomNumber = threshhold + 1e-6;
    bool isRadiating = mover->chooseDiagonalOrRadiating(
            sampler->getMovingBeads(),
            nslice, sampler->getSuperCell(), nslice);
    ASSERT_TRUE(isRadiating);
}

TEST_F(EMARateMoverTest, testDiagonalDecisionForIdenticalPaths) {
    positioner->setIdenticalPaths(separation);
    double delta2 = 3 * separation * separation;
    double diagonalProbability = 1.0;
    double radiatingProbability =
            pow(calculateRadiatingProbability(delta2), 2);
    double threshhold = diagonalProbability
            / (diagonalProbability + coefficient * radiatingProbability);
    mover->nextRandomNumber = threshhold - 1e-6;
    bool isRadiating = mover->chooseDiagonalOrRadiating(
            sampler->getMovingBeads(),
            nslice, sampler->getSuperCell(), nslice);
    ASSERT_FALSE(isRadiating);
}

TEST_F(EMARateMoverTest, testRadiatingDecisionForIdenticalPaths) {
    positioner->setIdenticalPaths(separation);
    double delta2 = 3 * separation * separation;
    double diagonalProbability = 1.0;
    double radiatingProbability =
            pow(calculateRadiatingProbability(delta2), 2);
    double threshhold = diagonalProbability
            / (diagonalProbability + coefficient * radiatingProbability);
    mover->nextRandomNumber = threshhold + 1e-6;
    bool isRadiating = mover->chooseDiagonalOrRadiating(
            sampler->getMovingBeads(),
            nslice, sampler->getSuperCell(), nslice);
    ASSERT_TRUE(isRadiating);
}
}
