#include <gtest/gtest.h>
#include "util/propagator/PotentialGrid.h"
#include <cmath>
#include <complex>

namespace {
typedef std::complex<double> Complex;

class PotentialGridTest: public testing::Test {
protected:

    void SetUp() {
        size = 20;
        deltaX = 0.1;
        x0 = -1.0;
        v = PotentialGridTest::myPotential;
        deltaTau = 0.1;
        grid = new PotentialGrid(size, deltaX, x0, v, deltaTau);
    }

    static double myPotential(double x) {
        return sin(x);
    }

    int size;
    double deltaX;
    double x0;
    double deltaTau;
    double (*v)(double);
    PotentialGrid* grid;
};

TEST_F(PotentialGridTest, TestValue) {
    double expect = exp(-deltaTau * v(-1.0));
    ASSERT_NEAR(expect, (*grid)(0), 1e-12);
}

}
