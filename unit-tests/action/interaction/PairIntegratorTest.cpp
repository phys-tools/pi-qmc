#include <gtest/gtest.h>
#include <cmath>
#include "config.h"

#include "action/interaction/PairIntegrator.h"
#include "action/interaction/InverseCosh2Potential.h"


namespace {

class PairIntegratorTest: public testing::Test {
protected:
  PairIntegratorTest()
    : v(1.0, 1.0), integrator(0)
  {
  }

  virtual void SetUp() {
    double tau = 0.1, mu = 1.0, dr = 0.1, tol = 1e-3, intRange = 5.0;
    int norder = 0, maxiter = 3, nsegment = 2;
    integrator = new PairIntegrator(tau, mu, dr, norder, maxiter,
      v,  tol, nsegment, intRange);
  }

  virtual void TearDown() {
    delete integrator;
  }

  InverseCosh2Potential v;
  PairIntegrator* integrator;
};

TEST_F(PairIntegratorTest, testCreate) {
  const PairIntegrator::Array& u = integrator->getU();
  ASSERT_EQ(1, u.size());
  EXPECT_DOUBLE_EQ(1.0, u(0));
}

}
