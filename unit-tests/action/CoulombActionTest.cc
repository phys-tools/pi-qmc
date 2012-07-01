#include <gtest/gtest.h>

#include "action/coulomb/Coulomb1DLinkAction.h"

namespace {

class CoulombActionTest: public ::testing::Test {
protected:

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }

};

TEST_F(CoulombActionTest, testValuesAtOrigin) {
    ASSERT_DOUBLE_EQ(0.0,
            Coulomb1DLinkAction::calculate1DValueAtOrigin(0.0));
    ASSERT_DOUBLE_EQ(0.017717130566979017,
            Coulomb1DLinkAction::calculate1DValueAtOrigin(0.01));
    ASSERT_DOUBLE_EQ(0.17650980432591801,
            Coulomb1DLinkAction::calculate1DValueAtOrigin(0.1));
    ASSERT_DOUBLE_EQ(-0.017731958122632733,
            Coulomb1DLinkAction::calculate1DValueAtOrigin(-0.01));
    ASSERT_DOUBLE_EQ(-0.17799263565884202,
            Coulomb1DLinkAction::calculate1DValueAtOrigin(-0.1));

}

}
