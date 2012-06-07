#include <gtest/gtest.h>

#include "action/PrimColloidalAction.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class PrimColloidalActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        species.count = npart;
        simInfo.setTau(0.1);
        B1 = 10;
        B2 = 20;
        V_lig = 1;
        V_cdte = 2;
        V_cdse = 3;
        ndim = 3;
    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species;
    SimulationInfo simInfo;
    double B1;
    double B2;
    double V_lig;
    double V_cdte;
    double V_cdse;
    int ndim;
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

TEST_F(PrimColloidalActionTest, getActionDifferenceForIdenticalPathsIsZero) {
    PrimColloidalAction action(B1, B2, V_lig, V_cdte, V_cdse, simInfo, ndim, species);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}

TEST_F(PrimColloidalActionTest, getActionDifferenceForOneMovedBead) {
    PrimColloidalAction action(B1, B2, V_lig, V_cdte, V_cdse, simInfo, ndim, species);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.2, deltaAction);
}

//TEST_F(PrimColloidalActionTest, getActionDifferenceForOneMovedBead2) {
//    PrimColloidalAction action(B1, B2, V_lig, V_cdte, V_cdse, simInfo, ndim, species);
//    setIdenticalPaths();
//    Beads<NDIM> *movingBeads = sampler->movingBeads;
//    Beads<NDIM>::Vec position(10.0, 8.0, 1.0);
//    (*movingBeads)(0, 32) = position;
//    double deltaAction = action.getActionDifference(*sampler, 0);
//    ASSERT_FLOAT_EQ(0.3, deltaAction);
//}

}
