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
#include <blitz/tinyvec-et.h>

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
extern "C" void DGESV_F77(const int *n, const int *nrhs, double *a, 
 const int *lda, int *ipiv, const double *b, const int *ldb, int *info );

double CaoBerneAction::u(double r, int idata) const {
  if (r != q) calcU(r);
  return uarray(0,idata);
} 

double CaoBerneAction::utau(double r, int idata) const {
  if (r != q) calcU(r);
  return uarray(1,idata);
} 

CaoBerneAction::CaoBerneAction(double mu, double radius, double tau, int norder)
  : radius(radius), mu(mu), tau(tau), norder(norder),
    ndata((norder+1)*(norder+2)/2), uarray(2,ndata),
    mat(ndata,ndata), ipiv(ndata), s2(ndata), z(ndata) {
}

void CaoBerneAction::calcU(double qnew) const {
  q = qnew;
#if NDIM==3
  double dr = sqrt(tau/mu);
  if (dr>0.2*q/(norder+1)) dr = 0.2*q/(norder+1);
  //Now find coefficients U(s^2)=sum c_ij s^2j z^2i
  int idata = 0;
  for (int iorder=0; iorder <= norder; ++iorder) {
    for (int ioff=0; ioff <= iorder; ++ioff) {
      double r1 = q + (iorder-ioff)*dr;
      double r2 = q - (iorder-ioff)*dr;
      double theta = ioff*dr/q;
      double costheta = cos(theta);
      double cos2theta = cos(2*theta);
      double sintheta = sin(theta);
      uarray(0,idata) = getAction(r1,r2,cos2theta,tau);
      uarray(1,idata) = (getAction(r1,r2,cos2theta,1.01*tau)
                        -getAction(r1,r2,cos2theta,0.99*tau))/(0.02*tau);
      Vec v1( r1*sintheta,0,r1*costheta);
      Vec v2(-r2*sintheta,0,r2*costheta);
      s2(idata) = dot((v1-v2),(v1-v2));
      z(idata) = r1-r2;
      ++idata;
    }
  }
  int jdata = 0;
  for (int jorder=0; jorder <= norder; ++jorder) {
    for (int joff=0; joff <= jorder; ++joff) {
      for (int i=0; i<ndata; ++i) {
        mat(jdata,i) = pow(s2(i),jorder)*pow(z(i),2*joff);
      }
      ++jdata;
    }
  }
  int info, nrhs=2;
  DGESV_F77(&ndata,&nrhs,mat.data(),&ndata,ipiv.data(),
            uarray.data(),&ndata,&info);
#else
  uarray=0;
#endif
}

double CaoBerneAction::getAction(double r1, double r2, double costheta, 
    double dtau) const {
  if (r1<radius || r2<radius) return 50;
  return -log(1.-radius*((r1+r2)-radius)/(r1*r2)
                *exp(-(r1*r2+radius*radius-radius*(r1+r2))*(1+costheta)
                     *mu/dtau));
}
