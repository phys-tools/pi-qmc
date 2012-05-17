#include <gtest/gtest.h>
#include <action/Action.h>
#include "action/SHOAction.h"

#include "sampler/test/MultiLevelSamplerFake.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "Beads.h"


namespace {

class SHOActionTest: public ::testing::Test {
public:
    typedef blitz::TinyVector<double,NDIM> Vec;
protected:

    virtual void SetUp() {
        sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
        species.count = npart;
        simInfo.setTau(0.1);
        tau =0.1;
        omega = 1;
        mass = 1;
        center = (0,0,0);

    }

    virtual void TearDown() {
        delete sampler;
    }

    MultiLevelSamplerFake *sampler;
    Species species;
    SimulationInfo simInfo;
    double tau;
    double mass;
    double omega;
    Vec center;
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

TEST_F(SHOActionTest, getActionDifferenceForIdenticalPathsIsZero) {
    SHOAction action(tau, omega, mass, NDIM, species, center);
    setIdenticalPaths();
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.0, deltaAction);
}
TEST_F(SHOActionTest, getActionDifferenceForOneMovedBead) {
    SHOAction action(tau, omega, mass, NDIM, species, center);
    setIdenticalPaths();
    Beads<NDIM> *movingBeads = sampler->movingBeads;
    Beads<NDIM>::Vec position(1.0, 2.0, 3.0);
    (*movingBeads)(0, 32) = position;
    double deltaAction = action.getActionDifference(*sampler, 0);
    ASSERT_FLOAT_EQ(0.46635586, deltaAction);
}

}
