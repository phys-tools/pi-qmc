#include <gtest/gtest.h>
#include "sampler/CollectiveSectionMover.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <blitz/tinyvec-et.h>

typedef blitz::TinyVector<double, NDIM> Vec;

bool vecEquals(const Vec& r1, const Vec& r2) {
    return dot(r1-r2,r1-r2) < 1e-10;
}




namespace {

class CollectiveSectionMoverTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        radius = 1.0;
        amplitude = Vec(0.4, 0.0, 0.0);
        mover = new CollectiveSectionMover(radius, amplitude);
    }

    virtual void TearDown() {
        delete mover;
    }

    void ASSERT_VEC_EQ(const Vec& v1, const Vec& v2) const {
        ASSERT_FLOAT_EQ(v1(0), v2(0));
        ASSERT_FLOAT_EQ(v1(1), v2(1));
        ASSERT_FLOAT_EQ(v1(2), v2(2));
    }

    double radius;
    Vec amplitude;
    CollectiveSectionMover *mover;
};

TEST_F(CollectiveSectionMoverTest, testMoveAtCenter) {
    Vec oldr(0.0, 0.0, 0.0);
    Vec newr = mover->calcShift(oldr);
    Vec expect = Vec(0.4, 0.0, 0.0);
    ASSERT_VEC_EQ(expect, newr);
}


TEST_F(CollectiveSectionMoverTest, testMoveAwayFromCenter) {
    Vec oldr(0.0, 0.5, 0.1);
    Vec newr = mover->calcShift(oldr);
    Vec expect = Vec(0.296, 0.5, 0.1);
    ASSERT_VEC_EQ(expect, newr);
}

TEST_F(CollectiveSectionMoverTest, testDoesNotMoveOutsideOfRadius) {
    Vec oldr(0.0, 1.5, 0.1);
    Vec newr = mover->calcShift(oldr);
    Vec expect = Vec(0.0, 1.5, 0.1);
    ASSERT_VEC_EQ(expect, newr);
}
}

