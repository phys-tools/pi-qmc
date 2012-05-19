#include <gtest/gtest.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "util/SuperCell.h"
#include "sampler/DoubleCollectiveSectionSampler.h"
#include "BeadFactory.h"
#include "sampler/DoubleSectionChooser.h"

typedef blitz::TinyVector<double, NDIM> Vec;
typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;


namespace {

class CollectiveSectionSamplerTest: public ::testing::Test {
protected:

    virtual void SetUp() {
        levelCount = 3;
        particleCount = 5;

        sampler = 0;
        paths = 0;
        action = 0;
        doubleAction = 0;
        beadFactory = 0;
        sectionChooser = 0;
    }

    virtual void TearDown() {
        delete sampler;
        delete beadFactory;
        delete sectionChooser;
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
    Paths *paths;
    Action *action;
    DoubleAction *doubleAction;
    int levelCount;
    BeadFactory *beadFactory;
    DoubleSectionChooser *sectionChooser;
    int particleCount;


    void createSampler() {
        Paths *paths = 0;
        Action *action = 0;
        DoubleAction *doubleAction = 0;
        const BeadFactory *beadFactory = new BeadFactory();
        DoubleSectionChooser *sectionChooser =
                new DoubleSectionChooser(levelCount,particleCount,*paths,
                        *action,*doubleAction,*beadFactory);
//        sampler = new DoubleCollectiveSectionSampler(particleCount,
//                paths, sectionChooser,
//                action, doubleAction, beadFactory);
        delete beadFactory;
    }
};

//TEST_F(CollectiveSectionSamplerTest, testCreate) {
//    createSampler();
//}

//TEST_F(CollectiveSectionSamplerTest, testSizeOfMovingBeads) {
//    createSampler();
//    Beads<NDIM> *beads = &sampler->getMovingBeads();
//    int sliceCount = beads->getNSlice();
//    int movingCount = beads->getNPart();
//    ASSERT_EQ(1 << levelCount + 1, sliceCount);
//    ASSERT_EQ(maximumMovingCount, movingCount);
//}
}

