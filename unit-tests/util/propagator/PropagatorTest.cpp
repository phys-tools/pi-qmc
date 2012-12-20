#include <gtest/gtest.h>
#include "util/propagator/Propagator.h"
#include "util/propagator/PropagatorGrid.h"
#include <cmath>

namespace {

class PropagatorTest: public testing::Test {
protected:
    void SetUp() {
        mass = 1.0;
        omega = 1.0;
    }

    void TearDown() {
    }

    double K0(double x1, double x2, double tau) {
        const double PI = 3.141592653589793;
        double amplitude = 0.0;
        int size = 256;
        double deltaX = 0.1;
        for (int n = -size / 2; n < size / 2; ++n) {
            double kn = 2.0 * PI * n / (size * deltaX);
            double weight = exp(-0.5 * tau * kn * kn / mass);
            amplitude += weight * cos(kn * (x1 - x2));
        }
        amplitude /= size;
        return amplitude;
    }


    double K(double x1, double x2, double tau) {
        const double PI = 3.141592653589793;
        double sinhwt = sinh(omega * tau);
        double coshwt = cosh(omega * tau);
        return sqrt(mass * omega / (2.0 * PI * sinhwt))
                * exp(-(mass*omega*(x1*x1 + x2*x2)*coshwt - 2* x1*x2)
                        / (2.0 * sinhwt));
    }

    double mass;
    double omega;
};

TEST_F(PropagatorTest, TestCreate) {
    Propagator propagator;
    double value = propagator.evaluate();
    double expect = K0(1.0, 1.0, 0.124235);
    ASSERT_NEAR(expect, value, 1e-12);
}

TEST_F(PropagatorTest, TestSetupGrid) {
    Propagator propagator;
    propagator.setupGrid();
    PropagatorGrid* grid = propagator.getGrid();
    ASSERT_NE(grid, (PropagatorGrid*)0);
}

}
