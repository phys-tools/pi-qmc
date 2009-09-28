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


#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int *m, const int *n, double *a, 
                           const int *lda, int *ipiv, int *info);

#define DGETRS_F77 F77_FUNC(dgetrs,DGETRS)
extern "C" void DGETRS_F77(const char *trans, const int *n, const int *nrhs,
                           const double *a, const int *lda, const int *ipiv,
                           double *b, const int *ldb, int *info);

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
extern "C" void DGESV_F77(const int *n, const int *nrhs,
                          const double *a, const int *lda, int *ipiv,
                          double *b, const int *ldb, int *info);


OptEwaldSum::OptEwaldSum(const SuperCell& cell, int npart,
  double rcut, double kcut, double khalo, int npoly)
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
  IArray ipiv(npoly);
  Array ftemp(npoly);
  a=0.; b=0.;
  for (unsigned int ik=0; ik<kval.size(); ++ik) {
    // Compute Fourier transform of basis functions and 1/r.
    double k = kval[ik];
#if NDIM==3
    double v = -4.*PI/(k*k);
    Complex eikr = exp(Complex(0.,k*rcut));
    Complex temp = (eikr-1.)/Complex(0.,k);
    v += 4.*PI*imag(temp)/k;
    double rn = rcut;
    for (int n=0; n<2*npoly-1; ++n) {
      temp = (rn*eikr - (n+1.)*temp)/Complex(0.,k);
      rn *= rcut;
      if (n%2 == 0) {
        ftemp(n/2) = 4.*PI*imag(temp)/k;
      }
    }
#else
    double v = -2.*PI/k;
    if (k*rcut > 30.) {
      v = 0; 
      ftemp =0.;
    } else {
      double factor = -0.25*k*k*rcut*rcut;
      /// handle n=-1 term
      double term = 2.*PI;
      double fn = term*rcut;
      int m = 1;
      while (fabs(term)>1e-100) {
        term *= factor/(m*m);
        fn += term*rcut/(2*m+1);
        ++m;
      }
      v += fn;
      // loop over even n terms
      double rnplus2 = rcut*rcut;
      for (int n=0; n<2*npoly; n+=2) {
        double term = 2.*PI;
        ftemp(n/2) = term*rnplus2/(n+2); 
        int m = 1;
        while (fabs(term)>1e-100) {
          term *= factor/(m*m);
          ftemp(n/2) += term*rnplus2/(2*m+n+2);
          ++m;
        }
        rnplus2 *= rcut*rcut;
      }
    }
#endif
    for (int i=0; i<npoly; ++i) {
      b(i) += degen[ik]*v*ftemp(i);
      for (int j=0; j<npoly; ++j) {
        a(i,j) += degen[ik]*ftemp(i)*ftemp(j);
      }
    }
  }
  // Lagrange multipliers for value (c) and derivative (d) at rcut.
  int ncts=2;
  Array c(npoly), d(npoly);
  for (int i=0; i<npoly; ++i) {
    c(i) = pow(rcut,2*i);
    d(i) = (2*i)*pow(rcut,2*i-1);
  }
  // Get ready to solve  linear equations.
  int info, one=1;
  char trans='N';
  DGETRF_F77(&npoly,&npoly,a.data(),&npoly,ipiv.data(),&info);
  // Solve with no constraint (Lagrange multipliers are zero).
  Array t0(npoly), tc(npoly), td(npoly), t(npoly);
  t0 = b;
  DGETRS_F77(&trans,&npoly,&one,a.data(),&npoly,ipiv.data(),
             t0.data(),&npoly,&info);
  double f0 = evalFRcut(t0,0), df0 = evalFRcut(t0,1);
  // Guess Lagrange multipliers.
  double lamc = 10*rcut;
  double lamd = -rcut*rcut;
  tc = b + lamc*c;
  DGETRS_F77(&trans,&npoly,&one,a.data(),&npoly,ipiv.data(),
             tc.data(),&npoly,&info);
  td = b + lamd*d;
  DGETRS_F77(&trans,&npoly,&one,a.data(),&npoly,ipiv.data(),
             td.data(),&npoly,&info);
  // Linear extropolation gives great estimate of correct Lagrange multipliers.
  Array2 mat(ncts,ncts);
  Array rhs(ncts);
  IArray ipv(ncts);
  mat(0,0) = (evalFRcut(tc,0)-f0)/lamc;
  mat(0,1) = (evalFRcut(td,0)-f0)/lamd;
  mat(1,0) = (evalFRcut(tc,1)-df0)/lamc;
  mat(1,1) = (evalFRcut(td,1)-df0)/lamd;
  rhs(0) =  1./rcut - f0;
  rhs(1) = -1./(rcut*rcut) - df0;
  DGESV_F77(&ncts,&one,mat.data(),&ncts,ipv.data(),rhs.data(),&ncts,&info);
  lamc = rhs(0); lamd = rhs(1);
  // Solve with correct Lagrange multipliers.
  t = b + lamc*c + lamd*d;
  DGETRS_F77(&trans,&npoly,&one,a.data(),&npoly,ipiv.data(),
             t.data(),&npoly,&info);
  coef = t;
#endif
  setLongRangeArray();
  evalSelfEnergy();
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
  double fpioverk = 4.*PI*imag(temp)/k;
  v -= fpioverk;
  double rn = rcut;
  for (int n=0; n<2*npoly-1; ++n) {
    temp = (rn*eikr - (n+1.)*temp)/Complex(0.,k);
    rn *= rcut;
    if (n%2 == 0) {
      v += coef(n/2) * fpioverk;
    }
  }
  return v;
#else
  double v = 2.*PI/k;
  if (k*rcut < 30.) {
    /// handle n=-1 term
    double factor = -0.25*k*k*rcut*rcut;
    double term = 2.*PI;
    double fn = term*rcut;
    int m = 1;
    while (fabs(term)>1e-100) {
      term *= factor/(m*m);
      fn += term*rcut/(2*m+1);
      ++m;
    }
    v -= fn;
    /// sum over even n terms
    double rnplus2 = rcut*rcut;
    for (int n=0; n<2*npoly; n+=2) {
      double term = 2.*PI;
      double fn = term*rnplus2/(n+2); 
      int m = 1;
      while (fabs(term)>1e-100) {
        term *= factor/(m*m);
        fn += term*rnplus2/(2*m+n+2);
        ++m;
      }
      v += coef(n/2)*fn;
      rnplus2 *= rcut*rcut;
    }
  }
  return v;
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
  double fn = 2.*PI*rcut;
  double v = -fn;
  fn *= rcut/2.;
  for (int n=0; n<2*npoly; n+=2) {
    v += coef(n/2)*fn;
    fn *= (rcut*rcut)*(2+n)/(4+n);
  }
  return v;
#endif
}

double OptEwaldSum::evalFRcut(Array& t, int nd) {
  // Calculate value or nd-order derivative of f(rcut) using coeffients t.
  double v=0;
  for (int i=0; i<npoly; ++i) {
    double x=pow(rcut,2*i-nd);
    for (int j=0; j<nd; ++j) {
      x *= 2*i-j;
    }
    v += t(i)*x;
  }
  return v;
}
