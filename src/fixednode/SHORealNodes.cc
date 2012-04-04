//$Id$
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
#include <cstdlib>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "SHORealNodes.h"
#include "util/PeriodicGaussian.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "util/SuperCell.h"

#define ZGETRF_F77 F77_FUNC(zgetrf,ZGETRF)
extern "C" void ZGETRF_F77(const int*, const int*, std::complex<double>*,
                           const int*, const int*, int*);
#define ZGETRI_F77 F77_FUNC(zgetri,ZGETRI)
extern "C" void ZGETRI_F77(const int*, std::complex<double>*, const int*, 
                     const int*, std::complex<double>*, const int*, int*);

SHORealNodes::SHORealNodes(const SimulationInfo &simInfo,
  const Species &species, const double omega, const double temperature,
  const double b, const int maxlevel)
  : tau(simInfo.getTau()),temperature(temperature),mass(species.mass),
    charge(species.charge), omega(omega), b(b), npart(species.count),
    ifirst(species.ifirst),
    matrix((int)(pow(2,maxlevel)+0.1)+1),
    gradmat1(npart,npart), gradmat2(npart,npart),
    gradmattau(npart,npart),  ipiv(npart), lwork(npart*npart), work(lwork),
    det((int)(pow(2,maxlevel)+0.1)+1),
    omegac(charge*b/(mass*c)), omega1(sqrt(omega*omega+0.25*omegac*omegac)),
    sinh1(sinh(0.5*omega1/temperature)), cosh1(cosh(0.5*omega1/temperature)),
    sinhc(sinh(0.5*omegac/temperature)), coshc(cosh(0.5*omegac/temperature)),
    a1(mass*omega1/(4.0*sinh1)*(cosh1-coshc)),
    a2(mass*omega1/(4.0*sinh1)*(cosh1+coshc)), b1(mass*omega1*sinhc/sinh1),
    a1tau(mass*omega1/(4.0*sinh1*sinh1)*
		    (-omega1+omega1*cosh1*coshc-omegac*sinhc*sinh1)),
    a2tau(mass*omega1/(4.0*sinh1*sinh1)*
		    (-omega1-omega1*cosh1*coshc+omegac*sinhc*sinh1)),
    b1tau(mass*omega1/(sinh1*sinh1)*
		    (omegac*coshc*sinh1-omega1*cosh1*sinhc)) {
  std::cout << "SHORealNodes with temperature=" << temperature
            << ", mass=" << mass << ", omega=" << omega
            << ", and B=" << b << std::endl;
  for (unsigned int i=0; i<matrix.size(); ++i) {
    matrix[i] = new Matrix(npart,npart,ColMajor());
  }
}

SHORealNodes::~SHORealNodes() {
  for (unsigned int i=0; i<matrix.size(); ++i) delete matrix[i];
}

NodeModel::DetWithFlag
SHORealNodes::evaluate(const VArray &r1, const VArray &r2, 
                          const int islice, bool scaleMagnitude) {
  DetWithFlag result; result.err=false;
  // First evaluate the inverse slater matrix.
  Matrix& mat(*matrix[islice]);
  for(int jpart=0; jpart<npart; ++jpart) {
    for(int ipart=0; ipart<npart; ++ipart) {
      Vec sum=Vec(r1(jpart+ifirst)+r2(ipart+ifirst));
      Vec diff=Vec(r1(jpart+ifirst)-r2(ipart+ifirst));
      mat(ipart,jpart)= exp(-a1*dot(sum,sum)-a2*dot(diff,diff));
      //std::complex<double> m= exp(-a1*dot(sum,sum)-a2*dot(diff,diff)
      //  -Complex(0,1)*b1*(r1(ipart+ifirst)[0]*r2(jpart+ifirst)[1]
      //                   -r1(ipart+ifirst)[1]*r2(jpart+ifirst)[0]));
      }
  }
  int info=0;
  ZGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
  if (info!=0) {
    result.err = true;
    std::cout << "BAD RETURN FROM ZGETRF!!!!" << std::endl;
  }
  det(islice) = 1;
  for (int i=0; i<npart; ++i) {
    det(islice)*= mat(i,i); 
    det(islice) *= (i+1==ipiv(i))?1:-1;
  }
  ZGETRI_F77(&npart,mat.data(),&npart,ipiv.data(),work.data(),&lwork,&info);
  if (info!=0) {
    result.err = true;
    std::cout << "BAD RETURN FROM ZGETRI!!!!" << std::endl;
  }
  result.det = real(det(islice));
  return result;
}

void SHORealNodes::evaluateDistance(const VArray &r1, const VArray &r2, 
                         const int islice, Array& d1, Array& d2) { 
  const Matrix& mat(*matrix[islice]);
  double d2inv=0;
  for(int jpart=0; jpart<npart; ++jpart) {
    for(int ipart=0; ipart<npart; ++ipart) {
      Vec sum=Vec(r1(jpart+ifirst)+r2(ipart+ifirst));
      Vec diff=Vec(r1(jpart+ifirst)-r2(ipart+ifirst));
      std::complex<double> m= exp(-a1*dot(sum,sum)-a2*dot(diff,diff)
        -Complex(0,1)*b1*(r1(jpart+ifirst)[0]*r2(ipart+ifirst)[1]
                         -r1(jpart+ifirst)[1]*r2(ipart+ifirst)[0]));
      //Compiler bug (Jan 23, 2012)
      //      gradmat1(ipart,jpart)=m*(-2*a1*sum-2*a2*diff);
      //      gradmat2(ipart,jpart)=m*(-2*a1*sum+2*a2*diff);
    }
  }

  // Calculate the nodal distance from jpart particles from r1.
  for(int jpart=0; jpart<npart; ++jpart) {
    CVec logGrad = Complex(0); 
    for(int ipart=0; ipart<npart; ++ipart) {
       logGrad += mat(jpart,ipart)*gradmat1(ipart,jpart);
    }
    Vec reallogGrad=0.0;
    for (int i=0; i<NDIM; ++i) {
       reallogGrad[i] = real( det(islice)*logGrad[i] ) /real(det(islice));
    }
    d2inv+=dot(reallogGrad, reallogGrad);
  }
  
  // Calculate the nodal distance from ipart particles from r2.
  for(int ipart=0; ipart<npart; ++ipart) {
    CVec logGrad = Complex(0); 
    for(int jpart=0; jpart<npart; ++jpart) {
       logGrad -= mat(jpart,ipart)*gradmat2(ipart,jpart);
    }
    Vec reallogGrad=0.0;
    for (int i=0; i<NDIM; ++i) {
       reallogGrad[i] = real( det(islice)*logGrad[i] ) /real(det(islice));
    }       
    d2inv+=dot(reallogGrad, reallogGrad);
  }
  //return sqrt(2*mass/(d2inv*tau));       
}

void SHORealNodes::evaluateGradLogDist(const VArray &r1, const VArray &r2,
    const int islice, VMatrix &gradd1, VMatrix &gradd2,
    const Array &d1, const Array &d2) {
}

const double SHORealNodes::c(137.0359895);
const double SHORealNodes::pi(acos(-1));
