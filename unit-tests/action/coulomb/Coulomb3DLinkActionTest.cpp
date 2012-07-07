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

}

