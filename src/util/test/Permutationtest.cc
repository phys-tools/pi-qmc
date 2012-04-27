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

TEST_F(PermutationTest, testAssignmentOperater) {
    Permutation permutation(5);
    permutation[3]=4;
    permutation[4]=3;
    ASSERT_EQ(4, permutation[3]);

}

TEST_F(PermutationTest, testResetPermutation) {
    Permutation permutation(5);
    permutation[3]=4;
    permutation[4]=3;
    permutation.reset();
    for (int index=0; index<5;++index) {
        ASSERT_EQ(index, permutation[index]);
    }
}

TEST_F(PermutationTest, testCopyPermutation) {
    Permutation permutation(5);
    permutation[3]=4;
    permutation[4]=3;
    Permutation apermutation(permutation);
    ASSERT_EQ(4, apermutation[3]);

}

TEST_F(PermutationTest, testReversePermutation) {
    Permutation permutation(5);
    permutation[3]=4;
    permutation[4]=2;
    permutation[2]=3;
    Permutation apermutation(permutation);
    ASSERT_EQ(4, permutation[3]);
    ASSERT_EQ(2, permutation[4]);
    ASSERT_EQ(3, permutation[2]);
    ASSERT_EQ(1, permutation[1]);
    ASSERT_EQ(0, permutation[0]);
    permutation.setToInverse(apermutation);
    ASSERT_EQ(2, permutation[3]);
    ASSERT_EQ(3, permutation[4]);
    ASSERT_EQ(4, permutation[2]);
    ASSERT_EQ(0, permutation[0]);
    ASSERT_EQ(1, permutation[1]);
}

TEST_F(PermutationTest, testPrependPermutation) {
    Permutation permutation(5);
    permutation[3]=4;
    permutation[4]=1;
    permutation[2]=3;
    permutation[1]=2;
    Permutation apermutation(permutation);
    permutation.setToInverse(apermutation);
    permutation.prepend(apermutation);
    ASSERT_EQ(3, permutation[3]);
    ASSERT_EQ(4, permutation[4]);
    ASSERT_EQ(2, permutation[2]);
    ASSERT_EQ(0, permutation[0]);
    ASSERT_EQ(1, permutation[1]);


}

TEST_F(PermutationTest, testAppendPermutation) {
    Permutation permutation(5);
    permutation[3]=4;
    permutation[4]=1;
    permutation[2]=3;
    permutation[1]=2;
    Permutation apermutation(permutation);
    permutation.setToInverse(apermutation);
    permutation.append(apermutation);
    ASSERT_EQ(3, permutation[3]);
    ASSERT_EQ(4, permutation[4]);
    ASSERT_EQ(2, permutation[2]);
    ASSERT_EQ(0, permutation[0]);
    ASSERT_EQ(1, permutation[1]);


}



}

