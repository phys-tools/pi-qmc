#include <gtest/gtest.h>

#include "util/AperiodicGaussian.h"





namespace {

class AperiodicGaussianTest: public ::testing::Test {
protected:
};
TEST_F(AperiodicGaussianTest, testAperiodicGaussianValue) {
   AperiodicGaussian apg(1.0, 10, 2);
   double value = apg(0);
   ASSERT_EQ(1, value);
}

TEST_F(AperiodicGaussianTest, testAperiodicGaussianGradValue) {
   AperiodicGaussian apg(100, 1, 2);
   double gradvalue = apg.grad(0);
   ASSERT_EQ(0, gradvalue);
}

TEST_F(AperiodicGaussianTest, testEquality) {
   AperiodicGaussian apg1(1.0, 10, 2);
   AperiodicGaussian apg2(100,1,2);
   double value1 = apg1(0);
   double gradvalue1 = apg1.grad(0);
   double value2 = apg2(0);
   double gradvalue2 = apg2.grad(0);
   ASSERT_EQ(value2, value1);
   ASSERT_EQ(gradvalue1, gradvalue2);
}
}
