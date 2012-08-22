#include <gtest/gtest.h>
#include <cmath>
#include "config.h"

#include "action/interaction/SHOInteraction.h"

#include "advancer/MultiLevelSamplerFake.h"
#include "base/Beads.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"

namespace {

class SHOInteractionTest: public testing::Test {
protected:
    typedef Beads<NDIM>::Vec Vec;

    virtual void SetUp() {
        species1 = createSpecies();
        species2 = createSpecies();

        mu = 1.0 / (1.0 / species1->mass + 1.0 / species2->mass);
        deltaTau = 0.1;
        simInfo.setTau(deltaTau);
        omega = 1.0;

        nslice = 5;
        nlevel = 2;
    }

    virtual void TearDown() {
        delete sampler;
    }

    Species* createSpecies() {
        Species *species = new Species();
        species->count = 1;
        simInfo.speciesList.push_back(species);
        simInfo.speciesIndex.push_back(species);
        return species;
    }

    MultiLevelSamplerFake *sampler;
    Species *species1;
    Species *species2;
    SimulationInfo simInfo;
    double omega;
    double mu;
    double deltaTau;
    static const int npart = 2;
    static const int nmoving = 2;
    static int nlevel;
    static int nslice;

    void setPaths(Vec newPos1, Vec newPos2, Vec oldPos1, Vec oldPos2) {
        Beads<NDIM> &sectionBeads = sampler->getSectionBeads();
        Beads<NDIM> &movingBeads = sampler->getMovingBeads();
        for (int i = 0; i < nslice; ++i) {
            double newWeight = 1.0 - abs(i - nslice / 2) * 1.0 / (nslice / 2);
            sectionBeads(0, i) = oldPos1;
            sectionBeads(1, i) = oldPos2;
            movingBeads(0, i) = newWeight * newPos1
                    + (1.0 - newWeight) * oldPos1;
            movingBeads(1, i) = newWeight * newPos2
                    + (1.0 - newWeight) * oldPos2;
        }
    }

    double calculateSHOAction(double x1, double x2) {
        double sinhwt = sinh(omega * deltaTau);
        double coshwt = cosh(omega * deltaTau);
        double expect = 0.5 * mu * omega
                * ((x1 * x1 + x2 * x2) * coshwt - 2.0 * x1 * x2) / sinhwt;
        expect -= 0.5 * mu * (x1 - x2) * (x1 - x2) / deltaTau;
        return expect;
    }
};

int SHOInteractionTest::nslice;
int SHOInteractionTest::nlevel;

TEST_F(SHOInteractionTest, getActionDifferenceForIdenticalPathsIsZero) {
    sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
    SHOInteraction action(simInfo, omega, species1, species2);
    setPaths(Vec(0.0), Vec(0.0), Vec(0.0), Vec(0.0));
    double deltaAction = action.getActionDifference(*sampler, 0);
    double expect = 0.0;
    ASSERT_DOUBLE_EQ(expect, deltaAction);
}

TEST_F(SHOInteractionTest, getActionDifferenceForShortPaths) {
    nslice = 3;
    nlevel = 1;
    sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
    SHOInteraction action(simInfo, omega, species1, species2);
    double dist = 1.0;
    setPaths(Vec(dist), Vec(-dist), Vec(0.0), Vec(0.0));
    double deltaAction = action.getActionDifference(*sampler, 0);

    double x1 = 0.0;
    double x2 = 2.0 * sqrt(NDIM) * dist;
    double expect = calculateSHOAction(x1, x2);
    expect *= 2;
    ASSERT_NEAR(expect, deltaAction, 1e-12);
}

TEST_F(SHOInteractionTest, getActionDifferenceForLongerPaths) {
    sampler = new MultiLevelSamplerFake(npart, nmoving, nslice);
    SHOInteraction action(simInfo, omega, species1, species2);
    double dist = 1.0;
    setPaths(Vec(dist), Vec(-dist), Vec(0.0), Vec(0.0));
    double deltaAction = action.getActionDifference(*sampler, 0);

    double x1 = 0.0;
    double x2 = sqrt(NDIM) * dist;
    double x3 = 2.0 * x2;
    double expect = calculateSHOAction(x1, x2) + calculateSHOAction(x2, x3);
    expect *= 2;
    ASSERT_NEAR(expect, deltaAction, 1e-12);
}

}
