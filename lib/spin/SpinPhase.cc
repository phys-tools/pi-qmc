//$Id: SpinPhase.cc,v 1.5 2007/10/03 12:53:56 jshumwa Exp $
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
#include <blitz/tinyvec-et.h>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include "SpinPhase.h"
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

SpinPhase::SpinPhase(const SimulationInfo &simInfo,
  const Species &species, const double t,  const double bx, const double by,
  const double bz, const double gmubs, const int maxlevel)
  : SpinPhaseModel(simInfo.getNPart()), tau(simInfo.getTau()), 
    charge(species.charge), temperature(t), bx(bx), by(by), bz(bz), 
    gmubs(gmubs), tanhx(tanh(0.5*gmubs*bx/t)), tanhy(tanh(0.5*gmubs*by/t)),
    tanhz(tanh(0.5*gmubs*bz/t)), npart(species.count),
    ifirst(species.ifirst),  matrix((int)(pow(2,maxlevel)+0.1)+1),
    gradmat1(npart,npart), gradmat2(npart,npart), ipiv(npart),
    lwork(npart*npart), work(lwork), f(1) /* f(1.0/(16.0*pi*pi)) */  {
  std::cout << "SpinPhase with t=" << temperature << ", b=" << bx << " "
	    << by << " " << bz << " " 
            << ", gmubs=" << gmubs << std::endl;
  for (unsigned int i=0; i<matrix.size(); ++i) {
    matrix[i] = new Matrix(npart,npart,ColMajor());
  }
}

SpinPhase::~SpinPhase() {
  for (unsigned int i=0; i<matrix.size(); ++i) delete matrix[i];
}

const double SpinPhase::pi(acos(-1));

void SpinPhase::evaluate(const VArray &r1Array, const VArray &r2Array, 
                         const SArray &s1Array, const SArray &s2Array, 
                         const int islice) {
  gradPhi1=0.0; sgradPhi1=0.0; vecPot1=0.0;
  gradPhi2=0.0; sgradPhi2=0.0; vecPot2=0.0;
  for (int ipart=0; ipart<npart; ++ipart) {
    SVec s1=s1Array(ipart+ifirst);
    SVec s2=s2Array(ipart+ifirst);
    Complex rho=0;
    C4Vec gradrho1, gradrho2;
    // Identity contribution to density matrix.
    rho += Complex( s1[0]*s2[0]+s1[1]*s2[1]+s1[2]*s2[2]+s1[3]*s2[3],
                   -s1[0]*s2[3]+s1[1]*s2[2]-s1[2]*s2[1]+s1[3]*s2[0]);
    gradrho1[0]=Complex(+s2[0],-s2[3]);
    gradrho1[1]=Complex(+s2[1],+s2[2]);
    gradrho1[2]=Complex(+s2[2],-s2[1]);
    gradrho1[3]=Complex(+s2[3],+s2[0]);
    gradrho2[0]=Complex(+s1[0],+s1[3]);
    gradrho2[1]=Complex(+s1[1],-s1[2]);
    gradrho2[2]=Complex(+s1[2],+s1[1]);
    gradrho2[3]=Complex(+s1[3],-s1[0]);
    // Z-direction contribution to density matrix.
    rho += Complex( s1[0]*s2[0]-s1[1]*s2[1]-s1[2]*s2[2]+s1[3]*s2[3],
                   -s1[0]*s2[3]-s1[1]*s2[2]+s1[2]*s2[1]+s1[3]*s2[0])*tanhz;
    gradrho1[0]+=Complex(+s2[0],-s2[3])*tanhz;
    gradrho1[1]+=Complex(-s2[1],-s2[2])*tanhz;
    gradrho1[2]+=Complex(-s2[2],+s2[1])*tanhz;
    gradrho1[3]+=Complex(+s2[3],+s2[0])*tanhz;
    gradrho2[0]+=Complex(+s1[0],+s1[3])*tanhz;
    gradrho2[1]+=Complex(-s1[1],+s1[2])*tanhz;
    gradrho2[2]+=Complex(-s1[2],-s1[1])*tanhz;
    gradrho2[3]+=Complex(+s1[3],-s1[0])*tanhz;

    // X-direction contribution to density matrix.
    rho += Complex( s1[0]*s2[0]-s1[1]*s2[1]-s1[2]*s2[2]+s1[3]*s2[3],
                   -s1[0]*s2[3]-s1[1]*s2[2]+s1[2]*s2[1]+s1[3]*s2[0])*tanhx;
    gradrho1[0]+=Complex(+s2[0],-s2[3])*tanhx;
    gradrho1[1]+=Complex(-s2[1],-s2[2])*tanhx;
    gradrho1[2]+=Complex(-s2[2],+s2[1])*tanhx;
    gradrho1[3]+=Complex(+s2[3],+s2[0])*tanhx;
    gradrho2[0]+=Complex(+s1[0],+s1[3])*tanhx;
    gradrho2[1]+=Complex(-s1[1],+s1[2])*tanhx;
    gradrho2[2]+=Complex(-s1[2],-s1[1])*tanhx;
    gradrho2[3]+=Complex(+s1[3],-s1[0])*tanhx;

    // Y-direction contribution to density matrix.
    rho += Complex( s1[0]*s2[0]-s1[1]*s2[1]-s1[2]*s2[2]+s1[3]*s2[3],
                   -s1[0]*s2[3]-s1[1]*s2[2]+s1[2]*s2[1]+s1[3]*s2[0])*tanhy;
    gradrho1[0]+=Complex(+s2[0],-s2[3])*tanhy;
    gradrho1[1]+=Complex(-s2[1],-s2[2])*tanhy;
    gradrho1[2]+=Complex(-s2[2],+s2[1])*tanhy;
    gradrho1[3]+=Complex(+s2[3],+s2[0])*tanhy;
    gradrho2[0]+=Complex(+s1[0],+s1[3])*tanhy;
    gradrho2[1]+=Complex(-s1[1],+s1[2])*tanhy;
    gradrho2[2]+=Complex(-s1[2],-s1[1])*tanhy;
    gradrho2[3]+=Complex(+s1[3],-s1[0])*tanhy;

    // Phase gradient is the imaginary part of the log gradient.
    for (int i=0; i<4; ++i) {
      sgradPhi1(ipart+ifirst)[i]=imag(gradrho1[i]/rho);
      sgradPhi2(ipart+ifirst)[i]=imag(gradrho2[i]/rho);
    }
  }
}

void SpinPhase::evalOrbit(const Vec& r1, const Vec& r2, 
    const SVec& s1, const SVec& s2, Complex &rho, Vec &gradrho1, Vec &gradrho2,
    C4Vec &sgradrho1, C4Vec &sgradrho2) {
  rho=0;
  // Identity contribution to density matrix.
  rho += Complex( s1[0]*s2[0]+s1[1]*s2[1]+s1[2]*s2[2]+s1[3]*s2[3],
                 -s1[0]*s2[3]+s1[1]*s2[2]-s1[2]*s2[1]+s1[3]*s2[0]);
  sgradrho1[0]=Complex(+s2[0],-s2[3]);
  sgradrho1[1]=Complex(+s2[1],+s2[2]);
  sgradrho1[2]=Complex(+s2[2],-s2[1]);
  sgradrho1[3]=Complex(+s2[3],+s2[0]);
  sgradrho2[0]=Complex(+s1[0],+s1[3]);
  sgradrho2[1]=Complex(+s1[1],-s1[2]);
  sgradrho2[2]=Complex(+s1[2],+s1[1]);
  sgradrho2[3]=Complex(+s1[3],-s1[0]);
  // Z-direction contribution to density matrix.
  rho += Complex( s1[0]*s2[0]-s1[1]*s2[1]-s1[2]*s2[2]+s1[3]*s2[3],
                 -s1[0]*s2[3]-s1[1]*s2[2]+s1[2]*s2[1]+s1[3]*s2[0])*tanhz;
  sgradrho1[0]+=Complex(+s2[0],-s2[3])*tanhz;
  sgradrho1[1]+=Complex(-s2[1],-s2[2])*tanhz;
  sgradrho1[2]+=Complex(-s2[2],+s2[1])*tanhz;
  sgradrho1[3]+=Complex(+s2[3],+s2[0])*tanhz;
  sgradrho2[0]+=Complex(+s1[0],+s1[3])*tanhz;
  sgradrho2[1]+=Complex(-s1[1],+s1[2])*tanhz;
  sgradrho2[2]+=Complex(-s1[2],-s1[1])*tanhz;
  sgradrho2[3]+=Complex(+s1[3],-s1[0])*tanhz;

  // X-direction contribution to density matrix.
  rho += Complex( s1[0]*s2[0]-s1[1]*s2[1]-s1[2]*s2[2]+s1[3]*s2[3],
                 -s1[0]*s2[3]-s1[1]*s2[2]+s1[2]*s2[1]+s1[3]*s2[0])*tanhx;
  sgradrho1[0]+=Complex(+s2[0],-s2[3])*tanhx;
  sgradrho1[1]+=Complex(-s2[1],-s2[2])*tanhx;
  sgradrho1[2]+=Complex(-s2[2],+s2[1])*tanhx;
  sgradrho1[3]+=Complex(+s2[3],+s2[0])*tanhx;
  sgradrho2[0]+=Complex(+s1[0],+s1[3])*tanhx;
  sgradrho2[1]+=Complex(-s1[1],+s1[2])*tanhx;
  sgradrho2[2]+=Complex(-s1[2],-s1[1])*tanhx;
  sgradrho2[3]+=Complex(+s1[3],-s1[0])*tanhx;
  // Y-direction contribution to density matrix.
  rho += Complex( s1[0]*s2[0]-s1[1]*s2[1]-s1[2]*s2[2]+s1[3]*s2[3],
                 -s1[0]*s2[3]-s1[1]*s2[2]+s1[2]*s2[1]+s1[3]*s2[0])*tanhy;
  sgradrho1[0]+=Complex(+s2[0],-s2[3])*tanhy;
  sgradrho1[1]+=Complex(-s2[1],-s2[2])*tanhy;
  sgradrho1[2]+=Complex(-s2[2],+s2[1])*tanhy;
  sgradrho1[3]+=Complex(+s2[3],+s2[0])*tanhy;
  sgradrho2[0]+=Complex(+s1[0],+s1[3])*tanhy;
  sgradrho2[1]+=Complex(-s1[1],+s1[2])*tanhy;
  sgradrho2[2]+=Complex(-s1[2],-s1[1])*tanhy;
  sgradrho2[3]+=Complex(+s1[3],-s1[0])*tanhy;
}
