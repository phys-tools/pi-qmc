#include <gtest/gtest.h>

#include "action/coulomb/CoulombLinkAction.h"

namespace {

class CoulombLinkActionTest: public testing::Test {
protected:
	CoulombLinkAction* action;
	double q1q2;
	double epsilon;
	double deltaTau;

	virtual void SetUp() {
		q1q2 = -1.0;
		epsilon = 12.0;
		deltaTau = 0.1;
		action = new CoulombLinkAction(q1q2, epsilon, deltaTau);
	}

	virtual void TearDown() {
		delete action;
	}

};

TEST_F(CoulombLinkActionTest, testValuesAtOrigin) {
	CoulombLinkAction::Vec delta1 = 0.0;
	CoulombLinkAction::Vec delta2 = 0.0;
	double value = action->getValue(delta1, delta2);
	ASSERT_NEAR(-6.61375337e-02, value, 1e-8);
}

}
