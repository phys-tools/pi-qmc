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
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include "SHOPhase.h"
#include "PeriodicGaussian.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "SuperCell.h"

#define ZGETRF_F77 F77_FUNC(zgetrf,ZGETRF)
extern "C" void ZGETRF_F77(const int*, const int*, std::complex<double>*,
                           const int*, const int*, int*);
#define ZGETRI_F77 F77_FUNC(zgetri,ZGETRI)
extern "C" void ZGETRI_F77(const int*, std::complex<double>*, const int*, 
                     const int*, std::complex<double>*, const int*, int*);

SHOPhase::SHOPhase(const SimulationInfo &simInfo,
  const Species &species, const double omega, const double temperature,
  const double b, const int maxlevel)
  : PhaseModel(simInfo.getNPart()),
    tau(simInfo.getTau()),temperature(temperature), mass(species.mass),
    charge(species.charge), omega(omega), b(b), npart(species.count),
    ifirst(species.ifirst),  matrix((int)(pow(2,maxlevel)+0.1)+1),
    gradmat1(npart,npart), gradmat2(npart,npart),
    ipiv(npart), lwork(npart*npart), work(lwork),
    omegac(charge*b/(2.0*mass*c)), omega1(sqrt(omega*omega+omegac*omegac)),
    sinh1(sinh(0.5*omega1/temperature)), cosh1(cosh(0.5*omega1/temperature)),
    sinhc(sinh(0.5*omegac/temperature)), coshc(cosh(0.5*omegac/temperature)),
    a1(mass*omega1/(4.0*sinh1)*(cosh1-coshc)),
    a2(mass*omega1/(4.0*sinh1)*(cosh1+coshc)), b1(mass*omega1*sinhc/sinh1) {
  std::cout << "SHOPhase with temperature=" << temperature
            << ", species=" << species 
	    << ", mass=" << mass << ", omega=" << omega
            << ", charge=" << charge << ", and B=" << b << std::endl;
  for (unsigned int i=0; i<matrix.size(); ++i) {
    matrix[i] = new Matrix(npart,npart,ColMajor());
  }
}

SHOPhase::~SHOPhase() {
  for (unsigned int i=0; i<matrix.size(); ++i) delete matrix[i];
}

const double SHOPhase::c(1); //(137.0359895);
const double SHOPhase::pi(acos(-1));

void SHOPhase::evaluate(const VArray &r1, const VArray &r2, 
                          const int islice) { 
  // First evaluate the inverse slater matrix.
  Matrix& mat(*matrix[islice]);
  for(int jpart=0; jpart<npart; ++jpart) {
    for(int ipart=0; ipart<npart; ++ipart) {
      
      Vec sum=r1(jpart+ifirst)+r2(ipart+ifirst);
      Vec diff=r1(jpart+ifirst)-r2(ipart+ifirst);
#if NDIM==3
      sum[2]=0, diff[2]=0;
#endif		
      std::complex<double> m=mat(ipart,jpart)
        =mass*omega1/(2*pi*sinh1)
	 *exp(-a1*(sum[0]*sum[0]+sum[1]*sum[1])
	      -a2*(diff[0]*diff[0]+diff[1]*diff[1])
          -Complex(0,1)*b1*(r1(jpart+ifirst)[0]*r2(ipart+ifirst)[1]
                           -r1(jpart+ifirst)[1]*r2(ipart+ifirst)[0]) );
           
      gradmat1(ipart,jpart)     = m*(-2*a1*sum-2*a2*diff);
      gradmat1(ipart,jpart)[0] -= m*Complex(0,1)*b1*r2(ipart+ifirst)[1];
      gradmat1(ipart,jpart)[1] += m*Complex(0,1)*b1*r2(ipart+ifirst)[0];
      gradmat2(ipart,jpart)     = m*(-2*a1*sum+2*a2*diff);
      gradmat2(ipart,jpart)[0] += m*Complex(0,1)*b1*r1(jpart+ifirst)[1];
      gradmat2(ipart,jpart)[1] -= m*Complex(0,1)*b1*r1(jpart+ifirst)[0];
#if NDIM==3
      gradmat1(ipart,jpart)[2]=0;
      gradmat2(ipart,jpart)[2]=0;
#endif
    }
  }
  // Calculate determinant and inverse.
  int info=0;
  ZGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
  if (info!=0) std::cout << "BAD RETURN FROM ZGETRF!!!!" << std::endl;
  std::complex<double> det=1;
  for (int i=0; i<npart; ++i) det *= mat(i,i) * ((i+1==ipiv(i))?1.:-1.);
  phi = arg(det);
  ZGETRI_F77(&npart,mat.data(),&npart,ipiv.data(),work.data(),&lwork,&info);
  if (info!=0) std::cout << "BAD RETURN FROM ZGETRI!!!!" << std::endl;
  gradPhi1=0.0; gradPhi2=0.0; vecPot1=0.0; vecPot2=0.0;
  // Calculate the gradient.
  for(int jpart=0; jpart<npart; ++jpart) {
    CVec loggrad1 = Complex(0); 
    CVec loggrad2 = Complex(0); 
    for(int ipart=0; ipart<npart; ++ipart) {
       loggrad1 += mat(jpart,ipart)*gradmat1(ipart,jpart);
       loggrad2 += mat(ipart,jpart)*gradmat2(jpart,ipart);
    }
    for (int idim=0; idim<NDIM; ++idim) {
      gradPhi1(jpart+ifirst)[idim]=imag(loggrad1[idim]);
      gradPhi2(jpart+ifirst)[idim]=imag(loggrad2[idim]);
    }
    vecPot1(jpart+ifirst)[0] = -0.5*b*r1(jpart+ifirst)[1]*charge;
    vecPot1(jpart+ifirst)[1] = +0.5*b*r1(jpart+ifirst)[0]*charge;
    vecPot2(jpart+ifirst)[0] = +0.5*b*r2(jpart+ifirst)[1]*charge;
    vecPot2(jpart+ifirst)[1] = -0.5*b*r2(jpart+ifirst)[0]*charge;
  }
}
