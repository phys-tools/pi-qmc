#include <gtest/gtest.h>

namespace {

class EstimatorParserTest: public ::testing::Test {
protected:

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }

};

TEST_F(EstimatorParserTest, testCreate) {
    ASSERT_FLOAT_EQ(0.0, 0.0);
}

}

