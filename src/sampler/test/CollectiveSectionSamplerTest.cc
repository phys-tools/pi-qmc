#include <gtest/gtest.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "SuperCell.h"
#include "sampler/DoubleCollectiveSectionSampler.h"

typedef blitz::TinyVector<double, NDIM> Vec;
typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;


namespace {

class CollectiveSectionSamplerTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        sampler = 0;
    }

    virtual void TearDown() {
        delete sampler;
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

    DoubleCollectiveSectionSampler* sampler;
};

TEST_F(CollectiveSectionSamplerTest, testCreate) {
    sampler = new DoubleCollectiveSectionSampler(0,0,0,0,0);
}

}

