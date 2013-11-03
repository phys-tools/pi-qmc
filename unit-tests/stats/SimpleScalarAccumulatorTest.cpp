#include <gtest/gtest.h>
#include <blitz/tinyvec-et.h>
#include "stats/SimpleScalarAccumulator.h"

namespace {

class SimpleScalarAccumulatorTest: public testing::Test {
protected:
    virtual void SetUp() {
        accumulator = new SimpleScalarAccumulator(0);
    }

    virtual void TearDown() {
        delete accumulator;
    }

    ScalarAccumulator* accumulator;
};

TEST_F(SimpleScalarAccumulatorTest, testCreate) {
}

}
