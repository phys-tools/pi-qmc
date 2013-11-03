#include <gtest/gtest.h>
#include <utility>

#include "action/coulomb/Coulomb3DLinkAction.h"
#include "action/coulomb/Coulomb1DLinkAction.h"

namespace {

class Coulomb3DLinkActionTest: public testing::Test,
        public testing::WithParamInterface<std::pair<double, double> > {
protected:
    double stau;
    Coulomb1DLinkAction* coulomb1D;
    Coulomb3DLinkAction* coulomb3D;

	virtual void SetUp() {
	    stau = 0.01;
	    coulomb1D = new Coulomb1DLinkAction(stau);
	    coulomb3D = new Coulomb3DLinkAction(*coulomb1D);
	}

	virtual void TearDown() {
	    delete coulomb3D;
	    delete coulomb1D;
	}

};

TEST_F(Coulomb3DLinkActionTest, testValueAtOrigin) {
	double reff = 0.0;
    double value = coulomb3D->calculateU0(reff);
	double expect = 0.017717130566979017;
	ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(Coulomb3DLinkActionTest, testU0) {
    double reff = 1.0;
    double value = coulomb3D->calculateU0(reff);
    double expect = 9.9843305484870431e-5;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(Coulomb3DLinkActionTest, testU1) {
    double reff = 1.0;
    double value = coulomb3D->calculateU1(reff);
    double expect = 8.5031825295702207e-6;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(Coulomb3DLinkActionTest, testU2) {
    double reff = 1.0;
    double value = coulomb3D->calculateU2(reff);
    double expect = 1.3429133944492796e-6;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(Coulomb3DLinkActionTest, testU3) {
    double reff = 1.0;
    double value = coulomb3D->calculateU3(reff);
    double expect = 1.637776980179888e-7;
    ASSERT_DOUBLE_EQ(expect, value);
}

TEST_F(Coulomb3DLinkActionTest, testU4) {
    double reff = 1.0;
    double value = coulomb3D->calculateU4(reff);
    double expect = 4.0165622775654946e-8;
    ASSERT_DOUBLE_EQ(expect, value);
}

}

