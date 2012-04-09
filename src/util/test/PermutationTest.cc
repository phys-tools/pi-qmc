#include <gtest/gtest.h>

#include "util/Permutation.h"




namespace {

class PermutationTest: public ::testing::Test {
protected:
};

TEST_F(PermutationTest, testInitializesToIdenty) {
    Permutation permutation(5);
    for (int index = 0; index < 5; ++index) {
        ASSERT_EQ(index, permutation[index]);
    }
}

}

