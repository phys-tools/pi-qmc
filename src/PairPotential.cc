// $Id: PairPotential.cc 50 2009-05-11 20:05:24Z john.shumwayjr $
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
