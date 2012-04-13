#include <gtest/gtest.h>

#include "util/Hungarian.h"


namespace {

class HungarianTest: public ::testing::Test {
protected:
};

TEST_F(HungarianTest, testMinusIdentityMatrix) {
    double matrix[3][3] = { {-1.0,  0.0,  0.0},
                            { 0.0, -1.0,  0.0},
                            { 0.0,  0.0, -1.0}};
    Hungarian hungarian(3);
    hungarian.solve(*matrix);
    ASSERT_EQ(-3.0, hungarian.getSum());
    double expect[3] = {0,1,2};
    for (int index = 0; index < 3; ++index) {
        ASSERT_EQ(expect[index], hungarian[index]);
    }
}


TEST_F(HungarianTest, testShiftedMinusIdentityMatrix) {
    double matrix[3][3] = { { 0.0, -1.0,  0.0},
                            { 0.0,  0.0, -1.0},
                            {-1.0,  0.0,  0.0}};
    Hungarian hungarian(3);
    hungarian.solve(*matrix);
    ASSERT_EQ(-3.0, hungarian.getSum());
    double expect[3] = {2,0,1};
    for (int index = 0; index < 3; ++index) {
        ASSERT_EQ(expect[index], hungarian[index]);
    }
}

TEST_F(HungarianTest, testFindSmallestSumIsTen) {
    double matrix[3][3] = { { 9.0, 5.0, 4.0},
                            { 3.0, 2.0, 2.0},
                            { 5.0, 4.0, 2.0}};
    Hungarian hungarian(3);
    hungarian.solve(*matrix);
    ASSERT_EQ(10.0, hungarian.getSum());
    double expect[3] = {1,0,2};
    for (int index = 0; index < 3; ++index) {
        ASSERT_EQ(expect[index], hungarian[index]);
    }
}
}
