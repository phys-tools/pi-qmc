#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define NDIM 3
#include "EMARateAction.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "MultiLevelSamplerInterface.h"
#include "Beads.h"
#include "SuperCell.h"
#include "SectionChooser.h"

class MultiLevelSamplerStub : public MultiLevelSamplerInterface {
public:
    MultiLevelSamplerStub() 
        :   sectionBeads(new Beads<NDIM>(npart, nslice)),
            movingBeads(new Beads<NDIM>(nmoving, nslice)),
            movingIndex(new IArray(nmoving)),
            superCell(blitz::TinyVector<double, NDIM>(1.0)) {
    }
    virtual ~MultiLevelSamplerStub() {
        delete movingBeads;
        delete sectionBeads;
        delete movingIndex;
    }
    virtual Beads<NDIM>& getSectionBeads() const {
        return *sectionBeads;
    }
    virtual Beads<NDIM>& getMovingBeads() const {
        return *movingBeads;
    }
    virtual IArray& getMovingIndex() const {
        return *movingIndex;
    }
    virtual SuperCell& getSuperCell() const {
        return superCell;
    }
    virtual int getFirstSliceIndex() const {
        return 0;
    }
    //virtual const SectionChooser& getSectionChooser() const;
    static const int npart = 2;
    static const int nmoving = 2;
    static const int nslice = 65;
    Beads<NDIM> *sectionBeads;
    Beads<NDIM> *movingBeads;
    IArray *movingIndex;
    mutable SuperCell superCell;
};


TEST(MathTest, TwoPlusTwoEqualsFour) {
    EXPECT_EQ(2 + 2, 4);
}

TEST(MathTest, TwoPlusThreeEqualsFive) {
    EXPECT_EQ(2 + 3, 5);
}


class EMARateActionTest : public testing::Test {
public:

    Species species1, species2;
    SimulationInfo simInfo;

    virtual void SetUp() {
    }

};

TEST_F(EMARateActionTest, CanConstruct) {
    EMARateAction action(simInfo, species1, species2, 1.0);
    MultiLevelSamplerStub sampler;
    action.testableGetActionDifference(sampler, 0); 
}
