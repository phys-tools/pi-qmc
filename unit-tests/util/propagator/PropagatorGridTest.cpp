#include <gtest/gtest.h>
#include "util/propagator/PropagatorGrid.h"
#include <cmath>

namespace {

class PropagatorGridTest: public ::testing::Test {
protected:

    void SetUp() {
        size = 32;
        grid = new PropagatorGrid(size);
    }

    int size;
    PropagatorGrid *grid;
};

TEST_F(PropagatorGridTest, TestInitialization) {
    int index0 = size / 2;
    grid->initialize(index0);
    ASSERT_NEAR(real((*grid)(index0)), 1.0, 1e-12);
    ASSERT_NEAR(real((*grid)(10)), 0.0, 1e-12);
}

}
