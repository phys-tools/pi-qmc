#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define NDIM 3
#include "EMARateAction.h"
#include "Species.h"
#include "SimulationInfo.h"

TEST(MathTest, TwoPlusTwoEqualsFour) {
    EXPECT_EQ(2 + 2, 4);
}

TEST(MathTest, TwoPlusThreeEqualsFive) {
    EXPECT_EQ(2 + 3, 5);
}


class EMARateActionTest : public testing::Test {
public:

  Species species1, species2;
  SimulationInfo simInfo;

  virtual void SetUp() {
  }

};

TEST_F(EMARateActionTest, CanConstruct) {
  EMARateAction(simInfo, species1, species2, 1.0);
}
