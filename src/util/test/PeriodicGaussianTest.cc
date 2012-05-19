#include <gtest/gtest.h>
#include "util/PeriodicGaussian.h"

namespace {

class PeriodicGaussianTest: public testing::Test
{
protected:
    void SetUp() {
        alpha = 2.5;
        length = 1.5;
        gridCount = 10000;
        pg = new PeriodicGaussian(alpha, length, gridCount);
    }

    void TearDown() {
        delete pg;
    }

    PeriodicGaussian* pg;
    double alpha;
    double length;
    int gridCount;
};

TEST_F(PeriodicGaussianTest, testNumberOfTerms)  {
    int nmax = PeriodicGaussian::numberOfTerms(alpha, length);
    ASSERT_EQ(5, nmax);
}


TEST_F(PeriodicGaussianTest, testNumberOfTermsForLargeBox)  {
    int nmax = PeriodicGaussian::numberOfTerms(alpha, 20.0);
    ASSERT_EQ(64, nmax);
}

TEST_F(PeriodicGaussianTest, testValueAtZero) {
    double value = (*pg)(0.0);
    double expect = 1.0072131266104112;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(PeriodicGaussianTest, testGradientAtZero) {
    double grad = pg->grad(0.0);
    double expect = 0.0;
    ASSERT_DOUBLE_EQ(expect, grad);
}

TEST_F(PeriodicGaussianTest, testSecondDerivativeAtZero) {
    double d2 = pg->d2(0.0);
    double expect = -4.6303272041148782;
    ASSERT_DOUBLE_EQ(expect, d2);
}

TEST_F(PeriodicGaussianTest, testValueNearEdge) {
    double value = (*pg)(0.4999 * length);
    double expect = 0.49012750361989998;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(PeriodicGaussianTest, testGradientNearEdge) {
    double grad = pg->grad(0.4999 * length);
    double expect = -0.00066637454689608252;
    ASSERT_DOUBLE_EQ(expect, grad);
}

TEST_F(PeriodicGaussianTest, testSecondDerivativeNearEdge) {
    double d2 = pg->d2(0.4999 * length);
    double expect = 4.442496431737875;
    ASSERT_DOUBLE_EQ(expect, d2);
}

}
