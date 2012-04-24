#include <gtest/gtest.h>

#include "util/PeriodicGaussian.h"





namespace {

class PeriodicGaussianTest: public ::testing::Test {
protected:
};
TEST_F(PeriodicGaussianTest, testPeriodicGaussianValue) {
   PeriodicGaussian pg1(100, 1, 1);
   double value = pg1(0);
   ASSERT_EQ(1, value);
}

TEST_F(PeriodicGaussianTest, testPeriodicGaussianGradValue) {
   PeriodicGaussian pg1(100, 1, 1);
   double gradvalue = pg1.grad(0);
   ASSERT_EQ(0, gradvalue);
}

TEST_F(PeriodicGaussianTest, testSecondDriv) {
   PeriodicGaussian pg1(100, 1, 1);
   PeriodicGaussian pg2(1, 10, 1);
   double secdrivalue1 = pg1.d2(0);
   double secdrivalue2 = pg2.d2(0);
   ASSERT_TRUE(secdrivalue1 - 100*secdrivalue2 < 1e-6);

}
}
