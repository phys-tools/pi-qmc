#include <gtest/gtest.h>
#include "sampler/CollectiveSectionMover.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

typedef blitz::TinyVector<double, NDIM> Vec;
typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;


bool vecEquals(const Vec& r1, const Vec& r2) {
    return dot(r1-r2,r1-r2) < 1e-10;
}

namespace {

class CollectiveSectionMoverTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        radius = 1.0;
        amplitude = Vec(0.4, 0.0, 0.0);
        center = Vec(3.0, 4.0, 5.0);
        sliceCount = 9;
        mover = new CollectiveSectionMover(radius, amplitude,
                center, sliceCount);
    }

    virtual void TearDown() {
        delete mover;
    }

    void ASSERT_VEC_EQ(const Vec& v1, const Vec& v2) const {
        ASSERT_FLOAT_EQ(v1(0), v2(0));
        ASSERT_FLOAT_EQ(v1(1), v2(1));
        ASSERT_FLOAT_EQ(v1(2), v2(2));
    }

    void ASSERT_MAT_EQ(const Mat& mat1, const Mat& mat2) const {
        ASSERT_FLOAT_EQ(mat1(0,0), mat2(0,0));
        ASSERT_FLOAT_EQ(mat1(0,1), mat2(0,1));
        ASSERT_FLOAT_EQ(mat1(0,2), mat2(0,2));
        ASSERT_FLOAT_EQ(mat1(1,0), mat2(1,0));
        ASSERT_FLOAT_EQ(mat1(1,1), mat2(1,1));
        ASSERT_FLOAT_EQ(mat1(1,2), mat2(1,2));
        ASSERT_FLOAT_EQ(mat1(2,0), mat2(2,0));
        ASSERT_FLOAT_EQ(mat1(2,1), mat2(2,1));
        ASSERT_FLOAT_EQ(mat1(2,2), mat2(2,2));
    }

    double radius;
    Vec amplitude;
    Vec center;
    int sliceCount;
    CollectiveSectionMover *mover;

    static Mat matFromData(double *data) {
        Mat matrix;
        for (int i=0; i<9; ++i) {
            matrix.data()[i] = data[i];
        }
        return matrix;
    }
};

TEST_F(CollectiveSectionMoverTest, testValueOfNslice) {
    ASSERT_EQ(9, mover->getSliceCount());
}

TEST_F(CollectiveSectionMoverTest, testMoveAtCenter) {
    Vec oldr = center;
    int sliceIndex = 4;
    Vec newr = mover->calcShift(oldr, sliceIndex);
    Vec expect = Vec(0.4, 0.0, 0.0) + center;
    ASSERT_VEC_EQ(expect, newr);
}


TEST_F(CollectiveSectionMoverTest, testMoveAwayFromCenter) {
    Vec oldr = Vec(0.0, 0.5, 0.1) + center;
    int sliceIndex = 4;
    Vec newr = mover->calcShift(oldr, sliceIndex);
    Vec expect = Vec(0.296, 0.5, 0.1) + center;
    ASSERT_VEC_EQ(expect, newr);
}

TEST_F(CollectiveSectionMoverTest, testDoesNotMoveOutsideOfRadius) {
    Vec oldr = Vec(0.0, 1.5, 0.1) + center;
    int sliceIndex = 4;
    Vec newr = mover->calcShift(oldr, sliceIndex);
    Vec expect = Vec(0.0, 1.5, 0.1) + center;
    ASSERT_VEC_EQ(expect, newr);
}

TEST_F(CollectiveSectionMoverTest, testDoesNotMoveAtFirstSlice) {
    Vec oldr = Vec(0.0, 0.5, 0.1) + center;
    int sliceIndex = 0;
    Vec newr = mover->calcShift(oldr, sliceIndex);
    Vec expect = Vec(0.0, 0.5, 0.1) + center;
    ASSERT_VEC_EQ(expect, newr);
}

TEST_F(CollectiveSectionMoverTest, testMoveAwayFromCenterSlice) {
    Vec oldr = Vec(0.0, 0.5, 0.1) + center;
    int sliceIndex = 2;
    Vec newr = mover->calcShift(oldr, sliceIndex);
    Vec expect = Vec(0.222, 0.5, 0.1) + center;
    ASSERT_VEC_EQ(expect, newr);
}

TEST_F(CollectiveSectionMoverTest, testJacobianMatrixAtCenter) {
    Vec oldr = center;
    int sliceIndex = 4;
    Mat jacobian = mover->calcJacobian(oldr, sliceIndex);
    double data[9] = {1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0};
    Mat expect = matFromData(data);
    ASSERT_MAT_EQ(expect, jacobian);
}

TEST_F(CollectiveSectionMoverTest, testJacobianAwayFromCenter) {
    Vec oldr = Vec(0.0, 0.5, 0.1) + center;
    int sliceIndex = 4;
    Mat jacobian = mover->calcJacobian(oldr, sliceIndex);
    double data[9] = {1.0, -0.4, -0.08,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0};
    Mat expect = matFromData(data);
    ASSERT_MAT_EQ(expect, jacobian);
}

TEST_F(CollectiveSectionMoverTest, testJacobianAwayFromCenterSlice) {
    Vec oldr = Vec(0.0, 0.5, 0.1) + center;
    int sliceIndex = 2;
    Mat jacobian = mover->calcJacobian(oldr, sliceIndex);
    double data[9] = {1.0, -0.3, -0.06,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0};
    Mat expect = matFromData(data);
    ASSERT_MAT_EQ(expect, jacobian);
}

TEST_F(CollectiveSectionMoverTest, testJacobianOutsideOfRadius) {
    Vec oldr = Vec(0.0, 1.5, 0.1) + center;
    int sliceIndex = 2;
    Mat jacobian = mover->calcJacobian(oldr, sliceIndex);
    double data[9] = {1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0};
    Mat expect = matFromData(data);
    ASSERT_MAT_EQ(expect, jacobian);
}

TEST_F(CollectiveSectionMoverTest, testReverseMove) {
    Vec oldr = Vec(0.1, 0.3, -0.2) + center;
    int sliceIndex = 3;
    Vec newr = mover->calcShift(oldr, sliceIndex);
    Vec backr = mover->calcInverseShift(newr, sliceIndex);
    ASSERT_VEC_EQ(backr, oldr);
}

}

