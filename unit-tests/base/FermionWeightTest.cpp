#include <gtest/gtest.h>

#include "base/FermionWeight.h"
#include "base/SerialPaths.h"
#include "base/BeadFactory.h"
#include "util/SuperCell.h"


namespace {

class FermionWeightTest: public ::testing::Test {
protected:
    void SetUp() {
        weight = new FermionWeight();
        cell = new SuperCell(SuperCell::Vec(1.0, 1.0, 1.0));
    }

    void TearDown() {
        delete weight;
        delete cell;
    }

    Paths* createPaths(int npart) {
        double deltaTau = 0.01;
        int nslice = 10;
        return new SerialPaths(npart, nslice, deltaTau, *cell, beadFactory);
    }

    void setPermutation(const Permutation& permutation, Paths *paths) {
        Beads<NDIM> *beads = beadFactory.getNewBeads(paths->getNPart(), 2);
        paths->putBeads(1, *beads, permutation);
        delete beads;
    }

    FermionWeight *weight;
    BeadFactory beadFactory;
    SuperCell *cell;
};


TEST_F(FermionWeightTest, testPartitionCount) {
    ASSERT_EQ(weight->getPartitionCount(), 2);
}


TEST_F(FermionWeightTest, testOneParticle) {
    Paths *paths = createPaths(1);
    weight->evaluate(paths);
    ASSERT_NEAR(weight->getValue(1), 1.0, 1e-14);
    delete paths;
}

TEST_F(FermionWeightTest, testAPairOfPermutingParticles) {
    Paths *paths = createPaths(2);
    Permutation p(2);
    p[0] = 1; p[1] = 0;
    setPermutation(p, paths);
    weight->evaluate(paths);
    ASSERT_NEAR(weight->getValue(1), -1.0, 1e-14);
    delete paths;
}

}
