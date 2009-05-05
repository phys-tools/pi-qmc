// $Id: OptEwaldSum.cc 38 2009-04-09 20:01:17Z john.shumwayjr $
/*  Copyright (C) 2009 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "OptEwaldSum.h"
#include <cmath>
#include <vector>
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"
#include "Paths.h"
#include "Beads.h"
#include "MultiLevelSampler.h"

OptEwaldSum::OptEwaldSum(const SuperCell& cell, int npart,
  double rcut, double kcut, double khalo, int npoly, int ncts)
  : EwaldSum(cell, npart, rcut, kcut), npoly(npoly), coef(npoly) {
#if NDIM==2 || NDIM==3
  // Fit coefficients to minimize chi2 error, constraining derivatives.
  // First find k values in the halo.
  double kcut2 = kcut*kcut;
  double khalo2 = khalo*khalo;
  std::vector<double> kval;
  std::vector<int> degen;
  kval.reserve(100); degen.reserve(100);
  Vec deltak,dk2;
  IVec cut, bigcut;
  for (int i=0; i<NDIM; ++i) {
    deltak[i] = 2*PI/cell[i];
    dk2[i] = deltak[i]*deltak[i];
  }
  for (int i=0; i<khalo/deltak[0]; ++i) {
    double kx2 = i*i*dk2[0];
    int idx = (i==0) ? 1 : 2;
    for (int j=0; j<sqrt(khalo2-kx2)/deltak[1]; ++j) {
#if NDIM==3
      double kx2ky2 = kx2 + j*j*dk2[1];
      int idxy = idx*((j==0) ? 1 : 2);
      for (int k=0; k<sqrt(khalo2-kx2ky2)/deltak[2]; ++k) {
        double k2 = kx2ky2 + k*k*dk2[2];
        int idegen = idxy*((k==0) ? 1 : 2);
#else
        double k2 = kx2 + j*j*dk2[1];
        int idegen = idx*((j==0) ? 1 : 2);
#endif
        if (k2>kcut2) {
          kval.push_back(sqrt(k2));
          degen.push_back(idegen);
        }
#if NDIM==3
      }
#endif
    }
  }
  // Now fit the coefficients.
  Array2 a(npoly,npoly);
  Array b(npoly);
  a=0.; b=0.;
  for (unsigned int i=0; i<kval.size(); ++i) {
  }
#endif
}

OptEwaldSum::~OptEwaldSum() {
}

double OptEwaldSum::evalFR(const double r) const {
  double v=0.;
  if (r>rcut) {
    v = 1./r;
  } else {
    double r2 = r*r;
    for (int i=npoly-1; i>=0; --i) {
      v *= r2;
      v += coef(i);
    }
  }
  return v;
}

double OptEwaldSum::evalFR0() const {
  return coef(0);
}

double OptEwaldSum::evalFK(const double k) const {
#if NDIM==3
  double v = 4.*PI/(k*k);
  Complex eikr = exp(Complex(0.,k*rcut));
  Complex temp = (eikr-1.)/Complex(0.,k);
  double fn = 4.*PI/k*imag(temp);
  v -= 4.*PI*imag(temp)/k;
  double rn = rcut;
  for (int n=0; n<2*npoly-1; ++n) {
    temp = (rn*eikr - (n+1.)*temp)/Complex(0.,k);
    rn *= rcut;
    if (n%2 == 0) {
      v += coef(n/2) * 4.*PI*imag(temp)/k;
    }
  }
  return v;
#else
  return 0.;
#endif
}

double OptEwaldSum::evalFK0() const {
#if NDIM==3
  double fn = 4.*PI*rcut*rcut/2.;
  double v = -fn;
  fn *= rcut*2./3.;
  for (int n=0; n<2*npoly; n+=2) {
    v += coef(n/2)*fn;
    fn *= (rcut*rcut)*(n+3.)/(n+5.);
  }
  return v;
#else
  return 0.;
#endif
}

