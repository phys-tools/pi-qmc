#include <gtest/gtest.h>

#include "util/math/VPolyFit.h"


namespace {

class VPolyFitTest: public ::testing::Test {
protected:
    double* xdata;
    double* ydata;
    int dataCount;
    int dimension;
};

TEST_F(VPolyFitTest, testFit) {
    dataCount = 3;
    dimension = 2;
    double tempx [] = {0.1, 0.2, 0.3};
    double tempy [] = {1.0, 3.0, 5.0,
                       1.0, 1.5, 2.0};
    xdata = tempx;
    ydata = tempy;
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_EQ(-1.0, solution[0]);
    ASSERT_EQ(0.5, solution[1]);
    delete vpolyfit;
}


}
