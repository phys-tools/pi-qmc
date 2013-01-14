#include <gtest/gtest.h>
#include "util/math/VPolyFit.h"
#include <cmath>

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
    double ydata [] = {1.0}; // y(x) = 1.0
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
    // y(x) = 5 x + 0.5
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
    // y(x) = 50x^2 - 10x + 1.5
    double ydata [] = {1.0, 1.5, 3.0};
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_DOUBLE_EQ(1.5, solution[0]);
    delete vpolyfit;
}

TEST_F(VPolyFitTest, TestManyPoints) {
    dataCount = 10;
    dimension = 1;
    double *xdata = new double [dataCount];
    for (int i = 0; i < dataCount; ++i) xdata[i] = 0.1 * (i + 0.5);
    double *ydata = new double [dataCount];
    for (int i = 0; i < dataCount; ++i) ydata[i] = 1.5 * sin(0.1*xdata[i]-9.3);
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* solution = vpolyfit->getSolution();
    ASSERT_NEAR(1.5 * sin(-9.3), solution[0], 1e-12);
    delete vpolyfit;
    delete [] xdata;
    delete [] ydata;
}


TEST_F(VPolyFitTest, testVectorFit) {
    dataCount = 3;
    dimension = 2;
    double xdata [] = {0.1, 0.2, 0.3};
    // y1(x) = 50x^2 - 10x + 1.5,
    // y2(x) = -50x^2 + 30x - 1.5
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

TEST_F(VPolyFitTest, TestConvergence) {
    dataCount = 3;
    dimension = 2;
    double xdata [] = {0.3, 0.2, 0.1};
    // y1(x) = 50x^2 - 10x + 1.5
    // linear approximation is y1(x) = 15x - 1.5
    // So last correction should be 3.0 from the last point.
    // y2(x) = 20x^2 - 5x
    // linear approximation is y1(x) = 10x - 1.8
    // So last correction should be 1.8 from the last point.
    double ydata [] = {3.0, 1.2,
                       1.5, 0.2,
                       1.0, -0.2};
    VPolyFit* vpolyfit = new VPolyFit(dataCount, dimension, xdata, ydata);
    vpolyfit->fit();
    const double* lastDelta = vpolyfit->getLastDelta();
    ASSERT_DOUBLE_EQ(3.0, lastDelta[0]);
    ASSERT_DOUBLE_EQ(1.8, lastDelta[1]);
    delete vpolyfit;
}

}
