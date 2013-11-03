#include <gtest/gtest.h>
#include "util/propagator/PropagatorGrid.h"
#include <cmath>
#include <complex>

namespace {
typedef std::complex<double> Complex;

class PropagatorGridTest: public ::testing::Test {
protected:

    void SetUp() {
        size = 32;
        deltaX = 0.2;
        grid = new PropagatorGrid(size, deltaX, 0.0);
    }

    int size;
    double deltaX;
    PropagatorGrid *grid;

    static const double PI;
    static const Complex I;

    double K0(double distance, double tau, double mass) {
        double amplitude = 0.0;
        for (int n = -size / 2; n < size / 2; ++n) {
            double kn = 2.0 * PI * n / (size * deltaX);
            double weight = exp(-0.5 * tau * kn * kn / mass);
            amplitude += weight * cos(kn * distance);
        }
        amplitude /= size;
        return amplitude;
    }

    /// Approximate grid propagator with continuum propagator.
    double approximateK0(double distance, double tau, double mass) {
        const double PI = 3.141592653589793;
        double delta2 = distance * distance;
        double prefactor = deltaX / sqrt(2 * PI * tau / mass);
        return prefactor * exp(-0.5 * mass * delta2 / tau);
    }

    static double potential(double x) {
        return 1.57;
    }
};

const double PropagatorGridTest::PI = 3.141592653589793;
const Complex PropagatorGridTest::I = Complex(0.0, 1.0);


TEST_F(PropagatorGridTest, TestSanityOfOurAnalyticK0Propagator) {
    double mass = 0.85;
    double offset = 0.2;
    double deltaTau = 0.14234;
    double expect = approximateK0(offset, deltaTau, mass);
    double value = K0(offset, deltaTau, mass);
    ASSERT_NEAR(expect, value, 1e-10);
}

TEST_F(PropagatorGridTest, TestDeltaK) {
    double deltaK = grid->getDeltaK();
    double expect = 2.0 * PI / (size * deltaX);
    ASSERT_NEAR(expect, deltaK, 1e-12);
}

TEST_F(PropagatorGridTest, TestInitialization) {
    int index0 = size / 2;
    grid->initialize(index0);
    ASSERT_NEAR(1.0, real((*grid)(index0)), 1e-12);
    ASSERT_NEAR(0.0, real((*grid)(10)), 1e-12);
}

TEST_F(PropagatorGridTest, TestTransformToKSpace) {
    int index0 = 13;
    double position = index0 * deltaX;
    grid->initialize(index0);
    grid->toKSpace();
    int indexk = 9;
    double k = indexk * grid->getDeltaK();
    Complex value = (*grid)(indexk);
    Complex expect = exp(-I * k * position) / sqrt(size);
    ASSERT_NEAR(real(expect), real(value), 1e-12);
    ASSERT_NEAR(imag(expect), imag(value), 1e-12);
}

TEST_F(PropagatorGridTest, TestTransformToKSpaceAndBack) {
    int index0 = 13;
    grid->initialize(index0);
    grid->toKSpace();
    grid->toRealSpace();
    ASSERT_NEAR(1.0, real((*grid)(index0)), 1e-12);
}

TEST_F(PropagatorGridTest, TestKineticEvolution) {
    int index0 = 13;
    int offset = 0;
    double mass = 1.0;
    double deltaTau = 0.14234;
    grid->setupKineticPropagator(mass, deltaTau);
    grid->initialize(index0);
    grid->toKSpace();
    grid->evolveTDeltaTau();
    grid->toRealSpace();
    double delta2 = (offset * deltaX) * (offset * deltaX);
    double expect = K0(offset * deltaX, deltaTau, mass);
    ASSERT_NEAR(expect, K0(0.0, deltaTau, mass), 1e-12);
    ASSERT_NEAR(expect, real((*grid)(index0 + offset)), 1e-12);
}


TEST_F(PropagatorGridTest, TestPotentialEvolution) {
    int index0 = 13;
    int offset = -2;
    double mass = 1.0;
    double deltaTau = 0.1;
    grid->setupPotentialPropagator(potential, deltaTau);
    grid->initialize(index0);
    grid->evolveVDeltaTau();
    double expect = exp(-deltaTau * potential(0.0));
    ASSERT_NEAR(expect, real((*grid)(index0)), 1e-12);
}

TEST_F(PropagatorGridTest, TestHalfPotentialEvolution) {
    int index0 = 13;
    int offset = -2;
    double mass = 1.0;
    double deltaTau = 0.1;
    grid->setupPotentialPropagator(potential, deltaTau);
    grid->initialize(index0);
    grid->evolveVHalfDeltaTau();
    double expect = exp(-0.5 * deltaTau * potential(0.0));
    ASSERT_NEAR(expect, real((*grid)(index0)), 1e-12);
}

}
