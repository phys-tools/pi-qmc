#include <gtest/gtest.h>
#include "util/propagator/Propagator1D.h"
#include "util/propagator/PropagatorGrid.h"
#include "PropagatorTestUtil.h"
#include <cmath>

namespace {

class Propagator1DTest: public testing::Test {
protected:
  void SetUp() {
    util.mass = mass = 1.0;
    util.omega = omega = 1.0;
    tau = 0.124235;
    x0 = 1.0;
  }

  double mass;
  double omega;
  double tau;
  double x0;

  PropagatorTestUtil util;
};

TEST_F(Propagator1DTest, TestDiagonalKineticEvolution) {
  Propagator1D propagator(mass, tau, x0);
  propagator.setPotential(Propagator1D::zeroPotential);
  double value = propagator.evaluate();
  double deltaX = propagator.getGridSpacing();
  double expect = util.approximateK0(x0, x0, tau, deltaX);
  ASSERT_NEAR(expect, value, 1e-12);
}

TEST_F(Propagator1DTest, TestDiagonalSHOEvolution) {
  Propagator1D propagator(mass, tau, x0);
  propagator.setPotential(Propagator1D::harmonicPotential);
  double value = propagator.evaluate();
  double deltaX = propagator.getGridSpacing();
  double expect = util.K(x0, x0, tau, deltaX);
  ASSERT_NEAR(expect, value, 1e-9);
}

}
