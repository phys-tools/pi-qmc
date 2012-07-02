#include <gtest/gtest.h>
#include <utility>

#include "action/coulomb/Coulomb1DLinkAction.h"

namespace {

class Coulomb1DLinkActionTest: public ::testing::Test,
		public ::testing::WithParamInterface<std::pair<double, double> > {
protected:

	virtual void SetUp() {
	}

	virtual void TearDown() {
	}

};

TEST_P(Coulomb1DLinkActionTest, testValueAtOrigin) {
	double stau = GetParam().first;
	double expect = GetParam().second;
	Coulomb1DLinkAction coulomb1D(stau);
	double value = coulomb1D.calculateValueAtOrigin();
	ASSERT_DOUBLE_EQ(expect, value)
			<< "(tau = " << stau << ")";
}

INSTANTIATE_TEST_CASE_P(TimestepScan, Coulomb1DLinkActionTest,
		::testing::Values(
				std::pair<double, double> (0.0, 0.0),
				std::pair<double, double> (0.01, 0.017717130566979017),
				std::pair<double, double> (0.1, 0.17650980432591801),
				std::pair<double, double> (-0.01, -0.017731958122632733),
				std::pair<double, double> (-0.1, -0.17799263565884202)));

}
