#include <gtest/gtest.h>
#include "util/math/VPolyFit.h"

namespace {

class VPolyFitTest: public ::testing::Test {
protected:
    int dataCount;
    int dimension;
};

TEST_F(VPolyFitTest, TestOnePoint) {
    dataCount = 1;
    dimension = 1;
    double xdata [] = {0.1};
    double ydata [] = {1.0};
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_EQ(1.0, solution[0]);
    delete vpolyfit;
}

TEST_F(VPolyFitTest, TestTwoPoints) {
    dataCount = 2;
    dimension = 1;
    double xdata [] = {0.1, 0.3};
    double ydata [] = {1.0, 2.0};
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_DOUBLE_EQ(0.5, solution[0]);
    delete vpolyfit;
}

TEST_F(VPolyFitTest, TestThreePoints) {
    dataCount = 3;
    dimension = 1;
    double xdata [] = {0.1, 0.2, 0.3};
    double ydata [] = {1.0, 1.5, 3.0};
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_DOUBLE_EQ(1.5, solution[0]);
    delete vpolyfit;
}

TEST_F(VPolyFitTest, testVectorFit) {
    dataCount = 3;
    dimension = 2;
    double xdata [] = {0.1, 0.2, 0.3};
    double ydata [] = {1.0, 1.0,
                       1.5, 2.5,
                       3.0, 3.0};
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_DOUBLE_EQ(1.5, solution[0]);
    ASSERT_DOUBLE_EQ(-1.5, solution[1]);
    delete vpolyfit;
}

}
