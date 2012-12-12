#include <gtest/gtest.h>
#include "util/propagator/PotentialGrid.h"
#include <cmath>
#include <complex>

namespace {
typedef std::complex<double> Complex;

class PotentialGridTest: public ::testing::Test {
protected:

    void SetUp() {
    }

};

TEST_F(PotentialGridTest, TestInitialization) {
    int size = 20;
    int deltaX = 0.1;
    double x0 = -1.0;
    double (*v)(double) = sin;
    PotentialGrid *grid = new PotentialGrid(size, deltaX, x0, v);
}

}
