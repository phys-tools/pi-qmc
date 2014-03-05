#include <gtest/gtest.h>
#include "util/propagator/Propagator1D.h"
#include "util/propagator/PropagatorGrid.h"
#include "util/propagator/PotentialGrid.h"
#include "action/interaction/AzizPotential.h"
#include "PropagatorTestUtil.h"
#include <cmath>

namespace {

class Propagator1DTest: public testing::Test {
protected:
  void SetUp() {
    mass = 1.0;
    omega = 1.0;
    tau = 0.124235;
    x0 = 1.0;
  }

  double mass;
  double omega;
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

TEST_F(Propagator1DTest, TestDiagonalKineticEvolution) {
  Propagator1D propagator(mass, tau, x0);
  propagator.setPotential(Propagator1D::zeroPotential);
  double value = propagator.evaluate();
  double deltaX = propagator.getGridSpacing();
  double expect = util.approximateK0(x0, x0, tau, deltaX, mass);
  ASSERT_NEAR(expect, value, 1e-12);
}

TEST_F(Propagator1DTest, TestDiagonalSHOEvolution) {
  Propagator1D propagator(mass, tau, x0);
  propagator.setPotential(Propagator1D::harmonicPotential);
  double value = propagator.evaluate();
  double deltaX = propagator.getGridSpacing();
  double expect = util.K(x0, x0, tau, deltaX, mass, omega);
  ASSERT_NEAR(expect, value, 1e-9);
}

TEST_F(Propagator1DTest, test_functor) {
  SHO_potential* sho = new SHO_potential();
  ASSERT_DOUBLE_EQ(3.125, (*sho)(2.5));
}

TEST_F(Propagator1DTest, test_SHO_potential_functor) {
  Propagator1D propagator(mass, tau, x0);
  SHO_potential* sho = new SHO_potential();
  ASSERT_DOUBLE_EQ(3.125, (*sho)(2.5));
  propagator.setPotential(sho);
  double value = propagator.evaluate();
  double deltaX = propagator.getGridSpacing();
  double expect = util.K(x0, x0, tau, deltaX, mass, omega);
  ASSERT_NEAR(expect, value, 1e-9);
}

TEST_F(Propagator1DTest, test_Aziz) {
  Aziz_potential* aziz = new Aziz_potential();
  const double K = 1/3.15774646e5, A=1/0.52917721092;
  ASSERT_DOUBLE_EQ(-10.896464664640893*K, (*aziz)(-3.0*A));
  Propagator1D propagator(2*1837.0, 1/(20.0*K), 3.0*A);
  propagator.setPotential(aziz);
  double value = propagator.evaluate();
  ASSERT_TRUE(value==value);
  ASSERT_NEAR(0.00098316676124997554, value, 1e-9);
}

}
