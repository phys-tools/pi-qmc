#include <gtest/gtest.h>
#include "util/propagator/Propagator1D.h"
#include "util/propagator/Propagator3D.h"
#include "util/propagator/PropagatorGrid.h"
#include "util/propagator/PotentialGrid.h"
#include "action/interaction/AzizPotential.h"
#include "PropagatorTestUtil.h"
#include <cmath>

namespace {

class Propagator3DTest: public testing::Test {
protected:
  void SetUp() {
  }

  double mass;
  double tau;
  double x0;

  PropagatorTestUtil util;

  struct SHO_potential : public PotentialGrid::functor
  {
    SHO_potential() : k(1.0) {}
    double operator()(double r) {return 0.5 * k * r * r;}
    double k;
  };

  struct Aziz_potential : public PotentialGrid::functor
  {
    AzizPotential aziz;
    double operator()(double r) {return aziz(fabs(r));}
  };
};

TEST_F(Propagator3DTest, test_Aziz) {
  const double K = 1/3.15774646e5, A=1/0.52917721092;
  x0 = 3.0 * A;
  tau = 1.0 / (20.0 * K);
  mass = 2.0*1822.888;
  Aziz_potential* aziz = new Aziz_potential();
  EXPECT_DOUBLE_EQ(-10.896464664640893*K, (*aziz)(x0));

  Propagator3D propagator0(mass, tau, x0);
  propagator0.setPotential(Propagator1D::zeroPotential);
  double value0 = propagator0.evaluate();
  double deltaX = propagator0.getGridSpacing();
  double expect0 = util.approximateK0(x0, x0, tau, deltaX, mass);
  EXPECT_NEAR(expect0, value0, 1e-12);

  Propagator3D propagator(mass, tau, x0);
  propagator.setPotential(aziz);
  double value = propagator.evaluate();
  EXPECT_TRUE(value==value); // Check for NAN
  EXPECT_NEAR(0.00097615769752803757, value, 1e-9);

  double u = -log(value / value0);
  EXPECT_NEAR(-0.018233141923916649, u, 1e-10);
}

}


