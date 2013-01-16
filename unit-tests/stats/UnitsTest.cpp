#include <gtest/gtest.h>
#include <blitz/tinyvec-et.h>
#include "stats/Units.h"

namespace {

class UnitsTest: public ::testing::Test {
    virtual void SetUp() {
        units = new Units("Ha", "a0");
    }

    virtual void TearDown() {
        delete units;
    }

protected:
    Units* units;
};

TEST_F(UnitsTest, testHartreeConversion) {
    double scale = units->getEnergyScaleIn("Ha");
    ASSERT_NEAR(1.0, scale, 1e-14);
}

TEST_F(UnitsTest, testBohrRadiusConversion) {
    double scale = units->getLengthScaleIn("a0");
    ASSERT_NEAR(1.0, scale, 1e-14);
}

TEST_F(UnitsTest, testElectronVoltConversion) {
    double scale = units->getEnergyScaleIn("eV");
    double expect = 1.0 / 27.2113845;
    ASSERT_NEAR(expect, scale, 1e-14);
}

TEST_F(UnitsTest, testNanometerConversion) {
    double scale = units->getLengthScaleIn("nm");
    double expect = 1.0 / 0.05291772086;
    ASSERT_NEAR(expect, scale, 1e-14);
}
}
