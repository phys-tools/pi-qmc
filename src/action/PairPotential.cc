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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PairPotential.h"

double PairPotential::getScatteringLength(double mu,
    double rmax, double dr) const {
  //Use Numerov method to integrate SE. a=r-psi(r)/psi'(r).
  double psi = dr, psim=0.;
  double f = 2*mu*operator()(dr), fm = 2*mu*operator()(0.);
  double r = dr;
  while (r<rmax) {
    r += dr;
    double fp = 2*mu*operator()(r);
    double psip = (2.*psi - psim + dr*dr*(fm*psim + 10*f*psi)/12.)
                 /(1.-dr*dr*fp/12.);
    psim=psi; fm=f;
    psi=psip; f=fp;
  }
  return (r-0.5*dr) - 0.5*dr*(psi+psim)/(psi-psim);
}

const double PairPotential::Aziz::A=1.9221529e5;
const double PairPotential::Aziz::alpha=10.73520708;
const double PairPotential::Aziz::beta=-1.89296514;
const double PairPotential::Aziz::C6=1.34920045;
const double PairPotential::Aziz::C8=0.41365922;
const double PairPotential::Aziz::C10=0.17078164;
const double PairPotential::Aziz::D=1.4135;
const double PairPotential::Aziz::eps=10.94/3.1577465e5;
const double PairPotential::Aziz::r0=0.2970/0.05291772086;
