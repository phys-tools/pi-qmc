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

TEST_F(Coulomb1DLinkActionTest, testU0) {
    double stau = 0.1;
    double reff = 1.0;
    Coulomb1DLinkAction coulomb1D(stau);
    double uOrigin = coulomb1D.calculateValueAtOrigin();
    double value = coulomb1D.calculateU0(uOrigin, reff);
    ASSERT_DOUBLE_EQ(0.010025994988629595, value);
}

TEST_F(Coulomb1DLinkActionTest, testU1) {
    double stau = 0.1;
    double reff = 1.0;
    Coulomb1DLinkAction coulomb1D(stau);
    double value = coulomb1D.calculateU1(reff);
    ASSERT_DOUBLE_EQ(0.00085593234703468702, value);
}

TEST_F(Coulomb1DLinkActionTest, testU2) {
    double stau = 0.1;
    double reff = 1.0;
    Coulomb1DLinkAction coulomb1D(stau);
    double value = coulomb1D.calculateU2(reff);
    ASSERT_DOUBLE_EQ(0.00012817165840892054, value);
}

TEST_F(Coulomb1DLinkActionTest, testU3) {
    double stau = 0.1;
    double reff = 1.0;
    Coulomb1DLinkAction coulomb1D(stau);
    double value = coulomb1D.calculateU3(reff);
    ASSERT_DOUBLE_EQ(2.3263720535382824e-5, value);
}

TEST_F(Coulomb1DLinkActionTest, testU4) {
    double stau = 0.1;
    double reff = 1.0;
    Coulomb1DLinkAction coulomb1D(stau);
    double value = coulomb1D.calculateU4(reff);
    ASSERT_DOUBLE_EQ(4.3433495124893814e-6, value);
}
