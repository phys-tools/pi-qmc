// $Id$
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
#ifndef __PairPotential_h_
#define __PairPotential_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
class Species;
#include <cstdlib>
#include <blitz/array.h>

/// Class for setting up empirical pair action between particles.
/// Has some subclasses that implement common pair potentials.
/// We may add methods for parsing formulas later.
/// @author John Shumway. 
class PairPotential {
public:
  virtual double operator()(double r) const=0;
  class InvCosh2;
  class LennardJones;
  class Aziz;
  double getScatteringLength(double mu, double rmax, double dr) const;
};

/// Empirical short range potential commonly used for trapped atoms.
class PairPotential::InvCosh2 : public PairPotential {
public: 
  InvCosh2(double v0, double kappa) : v0(v0), kappa(kappa) {}
  virtual double operator()(double r) const {
    double coshKappaR = cosh(kappa*r);
    return v0/(coshKappaR*coshKappaR);
  }
  double v0, kappa;
};


/// Lennard-Jones 6-12.
class PairPotential::LennardJones : public PairPotential {
public: 
  LennardJones(double epsilon, double sigma)
    : epsilon(epsilon), sigma(sigma) {}
  virtual double operator()(double r) const {
    double x = sigma/r;
    x *= x*x; x *= x; 
    return 4*epsilon*x*(x-1);
  }
  double epsilon, sigma;
};

/// Aziz potential (HFD-B2). 
/// R. A. Aziz et al., Mol. Phys. 77, 321 (1992).
class PairPotential::Aziz : public PairPotential {
  public: 
  virtual double operator()(double r) const {
    double x = r/r0;
    double v = (x<D) ? exp(-(D/x-1)*(D/x-1)) : 1.;
    double x2inv = 1./(x*x);
    v *= -x2inv*x2inv*x2inv*(C6+x2inv*(C8+x2inv*C10));
    v += A*exp(-x*(alpha-beta*x));
    v *= eps;
    return v;
  }
  static const double A;
  static const double alpha;
  static const double beta;
  static const double C6;
  static const double C8;
  static const double C10;
  static const double D;
  static const double eps;
  static const double r0;
};
#endif
