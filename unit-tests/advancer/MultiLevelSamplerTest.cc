#include <gtest/gtest.h>

namespace {

class MultiLevelSamplerTest: public ::testing::Test {
protected:

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }

};

TEST_F(MultiLevelSamplerTest, testCreate) {
    ASSERT_FLOAT_EQ(0.0, 0.0);
}

}

