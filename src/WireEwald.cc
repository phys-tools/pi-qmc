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
#include <cmath>
#include <iostream>
#include <fstream>
#include <blitz/tinyvec-et.h>
#include "SuperCell.h"
#include "WireEwald.h"

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);

using namespace blitz;

WireEwald::WireEwald(const SuperCell& cell, const double d, const double rmax,
                     const Vec fitWidth, const IVec gridSize)
  : cell(cell), d(d), rmax(rmax), gridSize(gridSize),
    dx(fitWidth[0] / double(gridSize[0]-1)), 
    dy(fitWidth[1] / double(gridSize[1]-1)),
#if NDIM==3
    dz(fitWidth[2] / double(gridSize[2]-1)),
    LL(1),
#endif
    NN(3), MM(1), 
#if NDIM==2
    v(gridSize[0]*gridSize[1]), 
    x(gridSize[0]*gridSize[1]), 
    y(gridSize[0]*gridSize[1]),
    polyCoeff((NN+1)*(MM+1))
#endif
#if NDIM==3
    v(gridSize[0]*gridSize[1]*gridSize[2]),
    x(gridSize[0]*gridSize[1]*gridSize[2]),
    y(gridSize[0]*gridSize[1]*gridSize[2]),
    z(gridSize[0]*gridSize[1]*gridSize[2]),
    polyCoeff((NN+1)*(MM+1)*(LL+1))
#endif
{
  double L = cell.a[0];
  int M = ceil(rmax / L);
  v = 0.; x = 0.; y = 0.;
  // Evaluate Coulomb interaction on the grid.
  for (int i=0; i<gridSize[0]; ++i) {
    for (int j=0; j<gridSize[1]; ++j) {
#if NDIM==2
      x(i*gridSize[1]+j) = dx * i;
      y(i*gridSize[1]+j) = dy * j;
      for (int k=-M; k<M+1; ++k) {
        if (k==0) continue;
        v(i*gridSize[1]+j) += 2./sqrt(pow(dx * i + k * L,2) + pow(dy * j,2))
                - 1./sqrt(pow(dx * i + k * L,2) + pow(dy * j + d,2))
                - 1./sqrt(pow(dx * i + k * L,2) + pow(dy * j - d,2));
      }
#endif
#if NDIM==3
      for (int k=0; k<gridSize[2]; ++k) {
        x(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) = dx * i;
        y(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) = dy * j;
        z(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) = dz * k;
        for (int kk=-M; kk<M+1; ++kk) {
          if (kk==0) continue;
          v(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) +=
                2./sqrt(pow(dx * i + k * L,2) + pow(dy * j,2) +pow(dz * k,2))
              - 1./sqrt(pow(dx * i + k * L,2) + pow(dy * j + d,2)+pow(dz,2))
              - 1./sqrt(pow(dx * i + k * L,2) + pow(dy * j - d,2)+pow(dz,2));
        }
      }
#endif
    }
  }
  // Fit v to find the coefficients for the polynomial.
  curveFit();
}

WireEwald::~WireEwald() {
}

double WireEwald::operator()(const Vec &r) {
#if NDIM==2
  return polynomial(NN,MM,r[0],r[1],polyCoeff);
#endif
#if NDIM==3
  return polynomial(NN,MM,LL,r[0],r[1],r[2],polyCoeff);
#endif
}

double WireEwald::getMaxDiff() {
#if NDIM==2
  Array1 f(gridSize[0]*gridSize[1]);
#endif
#if NDIM==3
  Array1 f(gridSize[0]*gridSize[1]*gridSize[2]);
#endif
  f = 0.;
  for (int i=0; i<gridSize[0]; ++i) 
    for (int j=0; j<gridSize[1]; ++j)
#if NDIM==2
      f(i*gridSize[1]+j) = polynomial(NN,MM,x(i*gridSize[1]+j),y(i*gridSize[1]+j),polyCoeff);
#endif
#if NDIM==3
      for (int k=0; k<gridSize[2]; ++k)
        f(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) = 
           polynomial(NN,MM,LL,x(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 y(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 z(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 polyCoeff);
#endif
  return max(abs(v-f)); 
}

void WireEwald::curveFit() {
  // Define parameters
  double tau = 1e-3;
  double epsilon1 = 1e-15, epsilon2 = 1e-15, epsilon3 = 1e-15;
  int KMAX = 100; 
#if NDIM==2
  Array1 p((NN+1)*(MM+1)), pnew((NN+1)*(MM+1));
#endif
#if NDIM==3
  Array1 p((NN+1)*(MM+1)*(LL+1)), pnew((NN+1)*(MM+1)*(LL+1));
#endif  
  p = 1.;
  // Evaluate the fitting function on the grid.
  Array1 f(v.shape()[0]), fnew(v.shape()[0]);
  f = 0.; fnew = 0.;
  for (int i=0; i<gridSize[0]; ++i) 
    for (int j=0; j<gridSize[1]; ++j)
#if NDIM==2
      f(i*gridSize[1]+j) = polynomial(NN,MM,x(i*gridSize[1]+j),y(i*gridSize[1]+j),p);
#endif
#if NDIM==3
      for (int k=0; k<gridSize[2]; ++k)
        f(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) = 
           polynomial(NN,MM,LL,x(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 y(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 z(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 p);
#endif
  // Calculate Jacobian on the grid.
  DArray J(v.shape()[0],p.shape()[0]);
  J = 0.;
  for (int i=0; i<v.shape()[0]; ++i) {
    for (int j=0; j<p.shape()[0]; ++j) {
#if NDIM==2
      int k = floor(j/(MM+1.));
      J(i,j) = pow(x(i),2*k) * pow(y(i),2*(j-k*(MM+1)));
#endif
#if NDIM==3
      int k = floor(j/((MM+1.)*(LL+1.)));
      int kk = floor(j/(LL+1.));
      J(i,j) = pow(x(i),2*k) * pow(y(i),2*(kk-k*(MM+1))) 
               * pow(z(i),2*(j-kk*(LL+1)));
#endif
    }
  }
  firstIndex firstindex;
  secondIndex secondindex;
  thirdIndex thirdindex;
  DArray A(J.shape()[1], J.shape()[1]);
  A = sum(J(thirdindex, firstindex) * J(thirdindex, secondindex), thirdindex);
  // Generate an identity matrix.
  DArray identity(A.shape());
  identity = 0.;
  for (int i=0; i<A.shape()[0]; ++i) 
    identity(i,i) = 1.;
  Array1 G(A.shape()[0]), epsilonP(v.shape());
  epsilonP = v - f;
  G = sum(J(secondindex, firstindex) * epsilonP(secondindex), secondindex);
  int count = 0;
  double nu = 2., mu = A(0,0);
  for (int i=1; i<A.shape()[0]; ++i)
    if (mu < A(i,i)) mu = A(i,i);
  mu *= tau;
  int dimMat = A.shape()[0];
  int lwork = dimMat*dimMat;
  Array1 work(lwork);
  IArray ipiv(dimMat);
  Array1 deltaP(dimMat);
  DArray mat(A.shape());
  double rhoe = 0.;
  bool stop = max(abs(G)) < epsilon1;
  while(!stop && (count < KMAX)) {
    ++count;
    do {
      mat = A + mu * identity;
      // LU decomposition
      int info = 0;
      DGETRF_F77(&dimMat, &dimMat, mat.data(), &dimMat, ipiv.data(), &info);
      if (info!=0) {
        std::cout<<"LU decomposition failed!"<<std::endl;
        std::exit(-1);
      }
      DGETRI_F77(&dimMat, mat.data(), &dimMat, ipiv.data(), work.data(),
                 &lwork, &info); 
      if (info!=0) {
        std::cout<<"Matrix inversion failed!"<<std::endl;
        exit(-1);
      }
      deltaP = sum(mat(firstindex, secondindex) * G(secondindex), secondindex);
      if (sum(deltaP*deltaP) <= epsilon2*epsilon2*sum(p*p)) 
        stop = true;
      else {
        pnew = p + deltaP;
        for (int i=0; i<gridSize[0]; ++i) 
          for (int j=0; j<gridSize[1]; ++j)
#if NDIM==2
            fnew(i*gridSize[1]+j) = polynomial(NN,MM,x(i*gridSize[1]+j),y(i*gridSize[1]+j),pnew);
#endif
#if NDIM==3
            for (int k=0; k<gridSize[2]; ++k)
              fnew(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k) = 
                 polynomial(NN,MM,LL,x(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                   y(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                   z(i*gridSize[1]*gridSize[2]+j*gridSize[2]+k),
                 p);
#endif
        rhoe = (sum(epsilonP*epsilonP) - 
                       sum((v - fnew)*(v-fnew))) / 
                       sum(deltaP * (mu* deltaP + G));
        if (rhoe > 0) {
          p = pnew;
          f = fnew;
          epsilonP = v - f;
          G = sum(J(secondindex, firstindex) * (epsilonP(secondindex)), secondindex);
          stop = ((max(abs(G)) < epsilon1) || (sum(epsilonP*epsilonP)) < epsilon3);
          mu *=(0.3333333333333333333 > 1-pow(2*rhoe-1, 3))?0.3333333333333333333:
               (1-pow(2*rhoe-1,3));
          nu = 2.;
        }
        else {
          mu *= nu;
          nu *= 2.;
        }
      }
    }while(!(rhoe > 0) && !stop);
  }
  polyCoeff = p;
  std::cout<<polyCoeff<<std::endl;
}

#if NDIM==2
inline double WireEwald::polynomial(int Nx, int My, double x, double y, const Array1& p) {
  double f = 0.0;
  for (int i=0; i<=Nx; ++i)
    for(int j=0; j<=My; ++j)
      f += p(i*(My+1)+j) * pow(x,2.*i) * pow(y,2.*j);
  return f;
}
#endif
#if NDIM==3
inline double WireEwald::polynomial(int Nx, int My, int Lz, double x, double y, 
                                    double z, const Array1& p) {
  double f = 0.0;
  for (int i=0; i<=Nx; ++i)
    for(int j=0; j<=My; ++j)
      for(int k=0; k<=Lz; ++k)
        f += p(i*(My+1)*(Lz+1)+j*(Lz+1)+k) 
             * pow(x,2.*i) * pow(y,2.*j) * pow(z,2.*k);
  return f;
}
#endif
