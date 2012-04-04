// $Id$
/*  Copyright (C) 2008-2009 John B. Shumway, Jr.

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
#include "EwaldSum.h"
#include <cmath>
#include "SimulationInfo.h"
#include "Species.h"
#include "util/SuperCell.h"
#include "Paths.h"
#include "Beads.h"
#include "sampler/MultiLevelSampler.h"
#include <blitz/tinyvec-et.h>
#ifdef _OPENMP
#include <omp.h>
#endif

EwaldSum::EwaldSum(const SuperCell& cell, const int npart,
                   const double rcut, const double kcut)
  : cell(cell), npart(npart), rcut(rcut), kcut(kcut),
    q(npart),
#if NDIM==3
    ikmax((int)(cell.a[0]*kcut/(2*PI)),(int)(cell.a[1]*kcut/(2*PI)),
          (int)(cell.a[2]*kcut/(2*PI))),
    deltak(2*PI*cell.b[0],2*PI*cell.b[1],2*PI*cell.b[2]),
    eikx(blitz::shape(ikmax[0]+1,npart)),
    eiky(blitz::shape(-ikmax[1],0), blitz::shape(2*ikmax[1]+1,npart)),
    eikz(blitz::shape(-ikmax[2],0), blitz::shape(2*ikmax[2]+1,npart)),
#endif
#if NDIM==2
    ikmax((int)(cell.a[0]*kcut/(2*PI)),(int)(cell.a[1]*kcut/(2*PI))),
    deltak(2*PI*cell.b[0],2*PI*cell.b[1]),
    eikx(blitz::shape(ikmax[0]+1,npart)),
    eiky(blitz::shape(-ikmax[1],0), blitz::shape(2*ikmax[1]+1,npart)),
#endif
    oneOver2V(0.5/product(cell.a))
{
#if (NDIM==3) || (NDIM==2)
  // Assume charges are all +1, can be reset if evalSelfEnergy is called again.
  q=1.0;
  // Calculate k-vectors.
  totk=0;
  for (int kx=0; kx<=ikmax[0]; ++kx) {
    double kx2=kx*kx*deltak[0]*deltak[0];
    for (int ky=((kx==0) ? 0 : -ikmax[1]); ky<=ikmax[1]; ++ky) {
#if NDIM==3
      double ky2=ky*ky*deltak[1]*deltak[1];
      for (int kz=((kx==0 && ky==0)? 0 : -ikmax[2]); kz<=ikmax[2]; ++kz) {
        double k2=kx2+ky2+kz*kz*deltak[2]*deltak[2];
#else
        double k2 = kx2 + ky*ky*deltak[1]*deltak[1];
#endif
        if (k2<kcut*kcut && k2!=0) ++totk;
#if NDIM==3
      }
#endif
    }
  }
  std::cout << "Ewald: totk=" << 2*totk << std::endl;
  vk.resize(totk);
  kvec.resize(totk);
  kvec2.resize(totk);
  // Calculate k-vectors (by symmetry, only loop over half of k-space).
  int ikvec=0;
  for (int kx=0; kx<=ikmax[0]; ++kx) {
    double kx2=kx*kx*deltak[0]*deltak[0];
    for (int ky=((kx==0) ? 0 : -ikmax[1]); ky<=ikmax[1]; ++ky) {
#if NDIM==3
      double ky2=ky*ky*deltak[1]*deltak[1];
      for (int kz=((kx==0 && ky==0)? 0 : -ikmax[2]); kz<=ikmax[2]; ++kz) {
        double k2=kx2+ky2+kz*kz*deltak[2]*deltak[2];
#else
        double k2 = kx2 + ky*ky*deltak[1]*deltak[1];
#endif
        if (k2<kcut*kcut && k2!=0) {
#if NDIM==3
          kvec(ikvec) = IVec(kx,ky,kz);
#else
          kvec(ikvec) = IVec(kx,ky);
#endif
          kvec2(ikvec) = k2;
          ++ikvec;
        }
#if NDIM==3
      }
#endif
    }
  }
  // Subclasses should be sure to call setLongRangeArray and
  // evalSelfEnergy in their constructors.
#endif
}

void EwaldSum::setLongRangeArray() {
#if (NDIM==3) || (NDIM==2)
  for (int ikvec=0; ikvec<totk; ++ikvec) {
    vk(ikvec)=evalFK(sqrt(kvec2(ikvec)));
  }
#endif
}

EwaldSum::~EwaldSum() {
}

double EwaldSum::evalLongRange(const VArray& r) const {

  double sum=0;
  // Set up the exponential tables
#if NDIM==3 || NDIM==2
  const Complex I(0,1);
    for (int ipart=0; ipart<npart; ++ipart) {
      eikx(0,ipart)=1.0;
      eiky(0,ipart)=1.0;
#if NDIM==3
      eikz(0,ipart)=1.0;
#endif
    } 
  if (ikmax[0]>0) {
    for (int ipart=0; ipart<npart; ++ipart) {
      eikx(1,ipart)=exp(I*deltak[0]*r(ipart)[0]);
    }
  } 
  if (ikmax[1]>0) {
    for (int ipart=0; ipart<npart; ++ipart) {
      eiky(1,ipart)=exp(I*deltak[1]*r(ipart)[1]);
      eiky(-1,ipart)=conj(eiky(1,ipart));
    }
  }
#if NDIM==3
  if (ikmax[2]>0) {
    for (int ipart=0; ipart<npart; ++ipart) {
      eikz(1,ipart)=exp(I*deltak[2]*r(ipart)[2]);
      eikz(-1,ipart)=conj(eikz(1,ipart));
    }
  } 
#endif
  for (int kx=2; kx<=ikmax[0]; ++kx) {
    for (int ipart=0; ipart<npart; ++ipart) {
      eikx(kx,ipart)=eikx(kx-1,ipart)*eikx(1,ipart);
    }
  }
  for (int ky=2; ky<=ikmax[1]; ++ky) {
    for (int ipart=0; ipart<npart; ++ipart) {
      eiky(ky, ipart)=eiky(ky-1, ipart)*eiky(1,ipart);
      eiky(-ky, ipart)=eiky(-ky+1, ipart)*eiky(-1,ipart);
    }
  }
#if NDIM==3
  for (int kz=2; kz<=ikmax[2]; ++kz) {
    for (int ipart=0; ipart<npart; ++ipart) {
      eikz(kz,ipart)=eikz(kz-1,ipart)*eikz(1,ipart);
      eikz(-kz,ipart)=eikz(-kz+1,ipart)*eikz(-1,ipart);
    }
  }
#endif
  // Sum long range action over the k-vectors.
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
#pragma omp for reduction(+:sum)
#endif
  for (int ikvec=0; ikvec<totk; ++ikvec) {
    Complex csum=0.;
    for (int jpart=0; jpart<npart; ++jpart) {
#if NDIM==3
      csum+=q(jpart)*eikx(kvec(ikvec)[0],jpart)*eiky(kvec(ikvec)[1],jpart)
	*eikz(kvec(ikvec)[2],jpart);  
#else
      csum+=q(jpart)*eikx(kvec(ikvec)[0],jpart)*eiky(kvec(ikvec)[1],jpart);
                 
#endif
    }

 
     sum+=2*vk(ikvec)*norm(csum);
  }
  } //end of omp parallel section
#endif
 
   return sum*oneOver2V + selfEnergy;
}

void EwaldSum::evalSelfEnergy() {
  double Q = sum(q);
  selfEnergy = -0.5*sum(q*q)*evalFR0() + oneOver2V*evalFK0()*Q*Q;
}

const double EwaldSum::PI=3.14159265358979;
