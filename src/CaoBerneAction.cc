// $Id$
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
#include "CaoBerneAction.h"
#include "PairPotential.h"
#include "Species.h"

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
extern "C" void DGESV_F77(const int *n, const int *nrhs, double *a, 
 const int *lda, int *ipiv, const double *b, const int *ldb, int *info );

double CaoBerneAction::u(double r, int idata) const {
  return 0.;
} 

double CaoBerneAction::utau(double r, int idata) const {
  return 0.;
} 

CaoBerneAction::CaoBerneAction(double mu, double tau, int norder)
  : radius(radius), mu(mu), tau(tau), norder(norder),
    ndata((norder+1)*(norder+2)/2) {
  //Now find coefficients U(s^2)=sum c_ij s^2j z^2i
/*  Array2 mat(ndata,ndata);
  IArray ipiv(ndata);
  u = log(g0)-log(g);
  int jdata = 0;
  for (int jorder=0; jorder<norder+1; ++jorder) {
    for (int joff=0; joff < ((NDIM>1)?jorder+1:1); ++joff) {
      for (int i=0; i<ndata; ++i) {
        mat(jdata,i) = pow(s2(i),jorder)*pow(z(i),2*joff);
      }
      ++jdata;
    }
  }
  int info, nrhs=1;
  DGESV_F77(&ndata,&nrhs,mat.data(),&ndata,ipiv.data(),u.data(),&ndata,&info); */
}

double CaoBerneAction::getAction(double r1, double r2, double costheta) const {
  return -log(1.-radius*((r1+r2)-radius)/(r1*r2)
                *exp(-(r1*r2+radius*radius*+radius*(r1+r2))*(1+costheta)
                     *(0.5*mu/tau)));
}
