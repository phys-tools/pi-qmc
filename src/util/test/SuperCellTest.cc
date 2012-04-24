#include <gtest/gtest.h>
#include <blitz/tinyvec-et.h>
#include "util/SuperCell.h"


namespace {

class SuperCellTest: public ::testing::Test {
protected:
};

TEST_F(SuperCellTest, testSuperCellRcut) {
    typedef blitz::TinyVector<double,NDIM> Vec;
    Vec a(1, 2, 4);
    SuperCell supercell(a);
    supercell.computeRecipricalVectors();
    ASSERT_EQ(1, supercell.b[0]);
    ASSERT_EQ(0.5, supercell.b[1]);
    ASSERT_EQ(0.25, supercell.b[2]);
    ASSERT_EQ(0.25, supercell.rcut2);
}

TEST_F(SuperCellTest, testSuperCellPbc) {
    typedef blitz::TinyVector<double,NDIM> Vec;
    Vec a(1, 2, 4);
    Vec v(0.1, 2, 2.5);
    SuperCell supercell(a);
    supercell.computeRecipricalVectors();
    supercell.pbc(v);
    ASSERT_EQ(0.1, v[0]);
    ASSERT_EQ(0, v[1]);
    ASSERT_EQ(-1.5, v[2]);
}

}
