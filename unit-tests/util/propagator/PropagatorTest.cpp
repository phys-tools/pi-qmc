#include <gtest/gtest.h>
#include "util/propagator/Propagator.h"
#include <cmath>

namespace {

class PropagatorTest: public ::testing::Test {
protected:

    double K(double x1, double x2, double tau) {
        double mass = 1.0;
        double omega = 1.0;
        const double PI = 3.141592653589793;
        double sinhwt = sinh(omega * tau);
        double coshwt = cosh(omega * tau);
        return sqrt(mass * omega / (2.0 * PI * sinhwt))
                * exp(-(mass*omega*(x1*x1 + x2*x2)*coshwt - 2* x1*x2)
                        / (2.0 * sinhwt));
    }

};

TEST_F(PropagatorTest, TestCreate) {
    Propagator propagator;
    double value = propagator.evaluate();
    double expect = K(1.0, 1.0, 1.0);
    ASSERT_NEAR(expect, value, 1e-12);
}

}
