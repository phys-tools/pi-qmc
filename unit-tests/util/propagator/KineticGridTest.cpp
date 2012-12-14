#include <gtest/gtest.h>
#include "util/propagator/KineticGrid.h"
#include <cmath>
#include <complex>

namespace {
typedef std::complex<double> Complex;

class KineticGridTest: public testing::Test {
protected:

    void SetUp() {
        gridCount = 16;
        deltaK = 0.1;
        mass = 2.1;
        deltaTau = 0.05;
        grid = new KineticGrid(gridCount, deltaK, mass, deltaTau);
    }

    void TearDown() {
        delete grid;
    }

    int gridCount;
    double deltaK;
    double mass;
    double deltaTau;
    KineticGrid* grid;
};


TEST_F(KineticGridTest, TestTheZeroKPoint) {
    double expect = 1.0;
    ASSERT_NEAR(expect, (*grid)(0), 1e-12);
}

TEST_F(KineticGridTest, TestAPositiveKPoint) {
    int index = 5;
    double k = index * deltaK;
    double expect = exp(-deltaTau * 0.5 * k * k / mass);
    ASSERT_NEAR(expect, (*grid)(index), 1e-12);
}

TEST_F(KineticGridTest, TestANegativeKPoint) {
    int index = 5;
    double k = - index * deltaK;
    double expect = exp(-deltaTau * 0.5 * k * k / mass);
    ASSERT_NEAR(expect, (*grid)(gridCount - index), 1e-12);
}

TEST_F(KineticGridTest, TestLargestKPoint) {
    int index = gridCount / 2;
    double k =  index * deltaK;
    double expect = exp(-deltaTau * 0.5 * k * k / mass);
    ASSERT_NEAR(expect, (*grid)(index), 1e-12);
}


}
