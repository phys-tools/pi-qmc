#include <gtest/gtest.h>
#include "util/PeriodicGaussian.h"

namespace {

class PeriodicGaussianTest: public testing::Test
{
protected:
    void SetUp() {
        alpha = 2.5;
        length = 1.5;
        pg = new PeriodicGaussian(alpha, length);
    }

    void TearDown() {
        delete pg;
    }

    PeriodicGaussian* pg;
    double alpha;
    double length;
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
    double value = pg->evaluate(0.0);
    double expect = 1.0072131266104112;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(PeriodicGaussianTest, testGradientAtZero) {
    pg->evaluate(0.0);
    double grad = pg->getGradient();
    double expect = 0.0;
    ASSERT_DOUBLE_EQ(expect, grad);
}

TEST_F(PeriodicGaussianTest, testSecondDerivativeAtZero) {
    pg->evaluate(0.0);
    double d2 = pg->getSecondDerivative();
    double expect = -4.6303272041148782;
    ASSERT_DOUBLE_EQ(expect, d2);
}

TEST_F(PeriodicGaussianTest, testValueNearEdge) {
    double value = pg->evaluate(0.4999 * length);
    double expect = 0.49012750361989998;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(PeriodicGaussianTest, testGradientNearEdge) {
    pg->evaluate(0.4999 * length);
    double grad = pg->getGradient();
    double expect = -0.00066637454689608252;
    ASSERT_DOUBLE_EQ(expect, grad);
}

TEST_F(PeriodicGaussianTest, testSecondDerivativeNearEdge) {
    pg->evaluate(0.4999 * length);
    double d2 = pg->getSecondDerivative();
    double expect = 4.442496431737875;
    ASSERT_DOUBLE_EQ(expect, d2);
}

TEST_F(PeriodicGaussianTest, testSymmetry) {
    double value = pg->evaluate(-0.3 * length);
    double expect = 0.66635569537953232;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(PeriodicGaussianTest, testAntisymmetryOfDerivative) {
    pg->evaluate(-0.3 * length);
    double grad = pg->getGradient();
    double expect = 1.0233852136355668;
    ASSERT_DOUBLE_EQ(expect, grad);
}

TEST_F(PeriodicGaussianTest, testSymmetryOfSecondDerivative) {
    pg->evaluate(-0.3 * length);
    double d2 = pg->getSecondDerivative();
    double expect = 1.4777706970511808;
    ASSERT_DOUBLE_EQ(expect, d2);
}

}
