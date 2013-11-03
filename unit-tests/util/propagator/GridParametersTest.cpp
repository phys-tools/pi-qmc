#include <gtest/gtest.h>
#include "util/propagator/GridParameters.h"
#include <cmath>
#include <complex>

namespace {
typedef std::complex<double> Complex;

class GridParametersTest: public testing::Test {
protected:

    void SetUp() {
        mass = 1.0;
        tau = 0.124235;
        x0 = 1.0;
        deltaX = 0.005;
        parameters = new GridParameters(mass, tau, x0, deltaX);

    }

    void TearDown() {
        delete parameters;
    }

    double mass;
    double tau;
    double x0;
    double deltaX;
    GridParameters *parameters;
};


TEST_F(GridParametersTest, TestThermalWidth) {
    double width = GridParameters::calculateThermalWidth(mass, tau);
    double expect = sqrt(tau/mass);
    ASSERT_NEAR(expect, width, 1e-12);
}

TEST_F(GridParametersTest, TestX0) {
    double xmin = parameters->getXMin();
    double deltaX = parameters->getDeltaX();
    int index0 = parameters->getIndex0();
    ASSERT_NEAR(x0, xmin + deltaX * index0, 1e-12);
}

TEST_F(GridParametersTest, TestXMinLowEnough) {
    double xmin = parameters->getXMin();
    double width = GridParameters::calculateThermalWidth(mass, tau);
    ASSERT_LT(xmin, x0 - 4.0 * width);
}


TEST_F(GridParametersTest, TestXMinNotTooLow) {
    double xmin = parameters->getXMin();
    double width = GridParameters::calculateThermalWidth(mass, tau);
    ASSERT_GT(xmin, x0 - 8.0 * width);
}

TEST_F(GridParametersTest, TestPowerTwoCeiling) {
    ASSERT_EQ(4096, GridParameters::powerTwoCeiling(3622.88));
}

}
