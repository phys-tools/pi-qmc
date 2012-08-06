#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "TradEwaldSum.h"
#include <cmath>
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "util/erf.h"
#include "util/SuperCell.h"

TradEwaldSum::TradEwaldSum(const SuperCell& cell, const int npart,
			   const double rcut, const double kcut, const double kappa)
  : EwaldSum(cell,npart,rcut,kcut), kappa(kappa) {
  setLongRangeArray();
  evalSelfEnergy();
}

TradEwaldSum::~TradEwaldSum() {
}

double TradEwaldSum::evalFR(const double r) const {
 return erf(kappa*r)/r;
}

double TradEwaldSum::evalFR0() const {
 return 2*kappa/sqrt(PI);
}

double TradEwaldSum::evalFK(const double k) const {
 return 4*PI*exp(-(k*k)/(4*kappa*kappa))/(k*k);
}

double TradEwaldSum::evalFK0() const {
 return -PI/(kappa*kappa);
}

