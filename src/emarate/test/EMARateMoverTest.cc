#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include "SimulationInfo.h"
#include "emarate/EMARateMover.h"
#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"

namespace {

class EMARateMoverTest: public testing::Test {
protected:


    virtual void SetUp() {
        separation = 5.0;
        setupSimulationInfo();
        coefficient = 10.0;
        mover = new EMARateMover(simInfo.tau, 1.0, 1.0, maxlevel, coefficient);
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
    }

    virtual void TearDown() {
        delete mover;
        delete sampler;
    }

    void setupSimulationInfo() {
        species1 = new Species("h", 1, 1.0, 0.0, 0, 0);
        species2 = new Species("e", 1, 1.0, 0.0, 0, 0);
        simInfo.nslice = 128;
        simInfo.npart = 2;
        simInfo.speciesList.resize(simInfo.npart);
        simInfo.speciesList[0] = species1;
        simInfo.speciesList[1] = species2;
        simInfo.speciesIndex.resize(simInfo.npart);
        simInfo.speciesIndex[0] = species1;
        simInfo.speciesIndex[1] = species2;
    }

    EMARateMover *mover;
    MultiLevelSamplerFake *sampler;

    Species *species1;
    Species *species2;
    double coefficient;
    static const int npart = 2;
    int nmoving;
    static const int maxlevel = 3;
    static const int nslice = 1 << maxlevel;
    SimulationInfo simInfo;
    double separation;
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

TEST_F(EMARateMoverTest, testRadiatingProbability) {

}

}
