#include <gtest/gtest.h>
#include <cmath>
#include "config.h"

#include "action/interaction/InverseCosh2Potential.h"


namespace {

class InverseCosh2PotentialTest: public testing::Test {
protected:

  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(InverseCosh2PotentialTest, testValues) {
  double v0 = 1.7, k = 3.0;
  InverseCosh2Potential v(v0, k);
  EXPECT_DOUBLE_EQ(v0, v(0.0));
  EXPECT_DOUBLE_EQ(v0 * pow(cosh(2.3 * k), -2), v(2.3));
}

}
