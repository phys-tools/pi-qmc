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
    double expect = K(1.0, 1.0, 1.0);
    ASSERT_NEAR(expect, value, 1e-12);
}

TEST_F(PropagatorTest, TestSetupGrid) {
    Propagator propagator;
    propagator.setupGrid();
    PropagatorGrid* grid = propagator.getGrid();
    ASSERT_NE(grid, (PropagatorGrid*)0);
}

}
