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
        deltaX = 0.1;
        grid = new PropagatorGrid(size, deltaX);
    }

    int size;
    double deltaX;
    PropagatorGrid *grid;

    static const double PI;
    static const Complex I;
};

const double PropagatorGridTest::PI = 3.141592653589793;
const Complex PropagatorGridTest::I = Complex(0.0, 1.0);


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

}
