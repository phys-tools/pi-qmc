#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include "SimulationInfo.h"
#include "emarate/EMARateMover.h"
#include "DeterministicEMARateMover.h"
#include "sampler/test/MultiLevelSamplerFake.h"
#include "EMARateTestBeadPositioner.h"
#include "Species.h"
#include <iostream>

namespace {

class EMARateMoverTest: public testing::Test {
protected:
    typedef DeterministicEMARateMover::Vec Vec;

    virtual void SetUp() {
        deltaTau = 0.4;
        nmoving = 2;
        separation = 1.0;
        setupSimulationInfo();
        coefficient = 10.0;
        mover = new DeterministicEMARateMover(deltaTau, 1.0, 1.0,
                maxlevel, coefficient);
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
        simInfo.tau = deltaTau;
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
    double deltaTau;
    static const int npart = 2;
    int nmoving;
    static const int maxlevel = 3;
    static const int nslice = (1 << maxlevel) + 1;
    SimulationInfo simInfo;
    double separation;
    static const double PI;

    double calculateRadiatingProbability(double delta2) {
        double reducedMass = 0.5;
        double tau = deltaTau * (nslice - 1) / 2;
        double sigma2 = tau / reducedMass;
        return exp(-0.5 * delta2 / sigma2);
    }

    double calculateDiagonalProbability(double delta2) {
        double mass = 1.0;
        double tau = deltaTau * (nslice - 1);
        double sigma2 = tau / mass;
        return exp(-0.5 * delta2 / sigma2);
    }

    bool vectorsAreClose(Vec v1, Vec v2, double threshold) {
        double delta2 = dot(v1-v2, v1-v2);
        return delta2 < threshold * threshold;
    }
};

const double EMARateMoverTest::PI = 3.141592653589793;

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
                    nslice, sampler->getSuperCell());
    double delta2 = 3 * separation * separation;
    double expect = pow(calculateRadiatingProbability(delta2), 2);
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testRadiatingProbabilityForRecombiningPaths) {
    positioner->setRecombiningPaths(separation);
    double probability =
            mover->calculateRadiatingProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell());
    double expect = 1.0;
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testDiagonalProbabilityForIdenticalPaths) {
    positioner->setIdenticalPaths(separation);
    double probability =
            mover->calculateDiagonalProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell());
    double expect = 1.0;
    ASSERT_DOUBLE_EQ(expect, probability);
}

TEST_F(EMARateMoverTest, testDiagonalProbabilityForRecombiningPaths) {
    positioner->setRecombiningPaths(separation);
    double probability =
            mover->calculateDiagonalProbability(sampler->getMovingBeads(),
                    nslice, sampler->getSuperCell());
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
            nslice, sampler->getSuperCell());
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
            sampler->getMovingBeads(), nslice, sampler->getSuperCell());
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
            sampler->getMovingBeads(), nslice, sampler->getSuperCell());
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
            sampler->getMovingBeads(), nslice, sampler->getSuperCell());
    ASSERT_TRUE(isRadiating);
}

TEST_F(EMARateMoverTest, testSampleRadiatingAtHighestLevel) {
    positioner->setIdenticalPaths(separation);
    mover->nextGaussianRandomNumber = Vec(0.3, 0.3, 0.3);
    mover->nextGaussianRandomNumber = Vec(0.0, 0.0, 0.0);

    Beads<NDIM> &movingBeads = sampler->getMovingBeads();
    mover->sampleRadiating(nslice/2, nslice, movingBeads);
    Vec radiatingPointBefore = movingBeads(0, nslice/2);
    Vec radiatingPointAfter = movingBeads(1, nslice/2);
    Vec holePosition = movingBeads(0,0);
    Vec electronPosition = movingBeads(1,0);
    double sigma2 = (nslice/2) * deltaTau / (1.0 + 1.0);
    Vec expect = 0.5 * (holePosition + electronPosition)
            + mover->nextGaussianRandomNumber * sqrt(sigma2);
    ASSERT_TRUE(vectorsAreClose(expect, radiatingPointBefore, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expect, radiatingPointAfter, 1e-12));
}

TEST_F(EMARateMoverTest, testSampleRadiatingAtLowestLevel) {
    positioner->setIdenticalPaths(separation);
    mover->nextGaussianRandomNumber = Vec(0.0, 0.0, 0.0);
    Beads<NDIM> &movingBeads = sampler->getMovingBeads();
    int nstride = nslice/2;
    while (nstride > 0) {
        mover->sampleRadiating(nstride, nslice, movingBeads);
        nstride >>= 1;
    }
    Vec firstPointAfterHole = movingBeads(0, 1);
    Vec firstPointBeforeHole = movingBeads(0, nslice - 2);
    Vec firstPointAfterElectron = movingBeads(1, 1);
    Vec firstPointBeforeElectron = movingBeads(1, nslice - 2);
    Vec thirdPointAfterHole = movingBeads(0, 3);
    Vec thirdPointBeforeHole = movingBeads(0, nslice - 4);
    Vec thirdPointAfterElectron = movingBeads(1, 3);
    Vec thirdPointBeforeElectron = movingBeads(1, nslice - 4);

    Vec holePosition = movingBeads(0,0);
    Vec electronPosition = movingBeads(1,0);
    double x = 1 / (nslice - 1.0);
    Vec expectFirstHole = (1.0 - x) * holePosition + x * electronPosition;
    Vec expectFirstElectron = x * holePosition + (1.0 - x) * electronPosition;
    x = 3 / (nslice - 1.0);
    Vec expectThirdHole = (1.0 - x) * holePosition + x * electronPosition;
    Vec expectThirdElectron = x * holePosition + (1.0 - x) * electronPosition;

    ASSERT_TRUE(vectorsAreClose(expectFirstHole, firstPointAfterHole, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectFirstHole, firstPointBeforeHole, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectFirstElectron, firstPointAfterElectron, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectFirstElectron, firstPointBeforeElectron, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectThirdHole, thirdPointAfterHole, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectThirdHole, thirdPointBeforeHole, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectThirdElectron, thirdPointAfterElectron, 1e-12));
    ASSERT_TRUE(vectorsAreClose(expectThirdElectron, thirdPointBeforeElectron, 1e-12));
}

TEST_F(EMARateMoverTest, testTransitionProbabilityForIdenticalPaths) {
    positioner->setIdenticalPaths(separation);
    Beads<NDIM> &sectionBeads = sampler->getSectionBeads();
    Beads<NDIM> &movingBeads = sampler->getMovingBeads();
    double probability = mover->calculateTransitionProbability(nslice/2,
            movingBeads, sectionBeads, nslice, sampler->getSuperCell());
    double expect = 0.0;
    ASSERT_DOUBLE_EQ(expect, probability);
}

//TEST_F(EMARateMoverTest, testTransitionProbabilityForRecombiningPaths) {
//    positioner->setRecombiningPaths(separation);
//    Beads<NDIM> &sectionBeads = sampler->getSectionBeads();
//    Beads<NDIM> &movingBeads = sampler->getMovingBeads();
//    double probability = mover->calculateTransitionProbability(nslice/2,
//            movingBeads, sectionBeads, nslice, sampler->getSuperCell());
//    double reverseProbability = mover->calculateDiagonalProbability(movingBeads,
//            nslice, sampler->getSuperCell());
//    double forwardProbability = mover->calculateRadiatingProbability(
//            movingBeads, nslice, sampler->getSuperCell());
//    double expect = log(forwardProbability / reverseProbability);
//    ASSERT_DOUBLE_EQ(expect, probability);
//}
}
