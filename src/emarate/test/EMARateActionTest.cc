#include <gtest/gtest.h>

#include "emarate/EMARateAction.h"
#include <iostream>
#include "util/SuperCell.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class EMARateActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        simInfo.nslice = 128;
        sampler->firstSliceIndex = simInfo.nslice - nslice/2;
        coefficient = 10.0;
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species1, species2;
    SimulationInfo simInfo;
    static const int npart=2;
    static const int nmoving=2;
    static const int nlevel=6;
    static const int nslice=8;
    double coefficient;

    void setIdenticalPaths() {
        Beads<NDIM>::Vec electronPosition(0.0);
        Beads<NDIM>::Vec holePosition(1.0);
        Beads<NDIM> &sectionBeads(sampler->getSectionBeads());
        Beads<NDIM> &movingBeads(sampler->getMovingBeads());
        for (int i = 0; i < nslice; ++i) {
            sectionBeads(0,i) = holePosition;
            sectionBeads(1,i) = electronPosition;
            movingBeads(0,i) = holePosition;
            movingBeads(1,i) = electronPosition;
        }
    }

    void setRecombiningPaths() {
        Beads<NDIM>::Vec beforePosition(0.0, 0.0, 0.0);
        Beads<NDIM>::Vec afterPosition(1.0, 1.0, 1.0);
        Beads<NDIM> &sectionBeads(sampler->getSectionBeads());
        Beads<NDIM> &movingBeads(sampler->getMovingBeads());
        for (int islice = 0; islice < nslice; ++islice) {
            sectionBeads(0,islice) = beforePosition;
            sectionBeads(1,islice) = afterPosition;
            if (islice < nslice/2) {
                movingBeads(0,islice) = beforePosition;
                movingBeads(1,islice) = beforePosition;
            } else if (islice == nslice/2) {
                movingBeads(0, islice) = beforePosition;
                movingBeads(1, islice) = afterPosition;
            } else {
                movingBeads(0,islice) = afterPosition;
                movingBeads(1,islice) = afterPosition;
            }
        }
    }
};

TEST_F(EMARateActionTest, testActionDifferenceForIdenticalPathsIsZero) {
    EMARateAction action(simInfo, species1, species2, coefficient);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

TEST_F(EMARateActionTest, testActionDifferenceForRecombiningPaths) {
    EMARateAction action(simInfo, species1, species2, coefficient);
    setRecombiningPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    double oldAction = -log(1 + coefficient * exp(-3.0));
    double newAction = -log(1 + coefficient * exp(+3.0));
    double expect = newAction - oldAction;
    ASSERT_FLOAT_EQ(expect, deltaAction);
}


}

