// $Id: EmpiricalInteraction.cc,v 1.1 2008/07/15 20:17:42 jshumwa Exp $
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

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
#include "EmpiricalInteraction.h"


EmpiricalInteraction::Cosh2Potential::Cosh2Potential(
  double v0, double kappa) 
  : v0(v0), kappa(kappa) {
}

double EmpiricalInteraction::Cosh2Potential::operator()(double r) const {
  double coshKappaR = cosh(kappa*r);
  return v0/(coshKappaR*coshKappaR);
} 

double EmpiricalInteraction::u(double r, int iorder) const {
  return (iorder==0) ? v(r)*tau : 0.;
} 

double EmpiricalInteraction::utau(double r, int iorder) const {
  return (iorder==0) ? v(r) : 0.;
} 


EmpiricalInteraction::EmpiricalInteraction(const Potential& v, const double tau)
  : v(v), tau(tau) {
}
