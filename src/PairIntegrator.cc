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
#include "PairIntegrator.h"
#include "PairPotential.h"
#include <cmath>
#include <blitz/tinyvec-et.h>
#include <iostream>
#include <fstream>

typedef std::complex<double> Complex;


extern "C" void dgetrf_(const int*, const int*, double*, const int*, const int*, int*);
extern "C" void dgetri_(const int*, double*, const int*, const int*, double*, const int*, int*);
extern "C" void dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, const double *b, const int *ldb, int *info );
extern "C" void dcopy_(const int *n, const double *x, const int *incx,
                                     const double *y, const int *incy);
extern "C" void dscal_(const int *n, const double *da, const double *dx,
  const int *incx);
extern "C" void daxpy_(const int *n, const double *da, const double *dx,
  const int *incx, const double *dy, const int *incy);
extern "C" void zcopy_(const int *n, const Complex *x, const int *incx,
                                     const Complex *y, const int *incy);
extern "C" void zdscal_(const int *n, const double *da, const Complex *dx,
  const int *incx);
extern "C" void zaxpy_(const int *n, const double *da, const Complex *dx,
  const int *incx, const Complex *dy, const int *incy);

PairIntegrator::PairIntegrator(double tau, double mu, double dr,
  int norder, int maxiter, const PairPotential &pot, double tol,
  int nsegment, double intRange) 
  : tau(tau), mu(mu), dr(dr), norder(norder), maxiter(maxiter),
    ndata(NDIM>1 ? (norder+1)*(norder+2)/2 : norder+1),
    nstep(maxiter), nstepInv(maxiter),
    g(ndata), g0(ndata), s2(ndata), z(ndata), u(ndata), pot(pot),
    tol(tol), nsegment(nsegment) {
  // Figure out the number of grid points.
  // (power of two larger than intRange*sigma).
  double sigma = sqrt(tau/mu);
  ngrid = 2*(int)exp2(ceil(log2((intRange*sigma/dr))));
  if (ngrid < 16) {
   ngrid=16;
   dr = 2*intRange*sigma/ngrid;
  }
  ngridN = pow(ngrid,NDIM);
  // Allocate the grids.
  expVtau.resize(ngridN);
  expHalfVtau.resize(ngridN);
  expTtau.resize(ngridN);
  expFullTtau.resize(ngridN);
  IVec size;
  for (int i=0; i<NDIM; ++i) size(i) = ngrid;
  vgrid.resize(size); vgrid=0.;
  tgrid.resize(size); tgrid=0.;
  IVecN1 size1;
  IVecN2 size2;
  size2(0) = maxiter;
  size1(0) = size2(1) = ndata;
  for (int i=0; i<NDIM; ++i) size1(i+1) = size2(i+2) = ngrid;
  psi.resize(size2);
  psi0.resize(size1);
  worka.resize(ndata,ngridN);
  workc.resize(maxiter,ndata,ngridN);
  workd.resize(maxiter,ndata,ngridN);
  best.resize(ndata,ngridN);
  err.resize(ndata,ngridN);
  // Set the kinetic energy.
  double dk=2*PI/(dr*ngrid);
  double over2mu = 0.5/mu;
  for (int i=0; i<ngrid; ++i) {
    int ishift = ((i+ngrid/2) % ngrid) - ngrid/2;
    double kx2 = ishift*ishift*dk*dk; 
#if NDIM>1
    for (int j=0; j<ngrid; ++j) {
      int jshift = ((j+ngrid/2) % ngrid) - ngrid/2;
      double ky2 = jshift*jshift*dk*dk; 
#if NDIM>2
      for (int k=0; k<ngrid; ++k) {
        int kshift = ((k+ngrid/2) % ngrid) - ngrid/2;
        double kz2 = kshift*kshift*dk*dk; 
        tgrid(i,j,k) = (kx2+ky2+kz2)*over2mu;
      }
#else
      tgrid(i,j) = (kx2+ky2)*over2mu;
#endif
    }
#else
    tgrid(i) = kx2*over2mu;
#endif
  }
  nstep(0)=2; nstep(1)=4; nstep(2)=6;
  for (int i=3; i<maxiter; ++i) nstep(i) = 2*nstep(i-2);
  nstepInv = 1./nstep;
  // Set up FFT's.
  fftw_complex *ptr = (fftw_complex*)psi.data();
  int istride = 1;
  IVec ngridArray = ngrid;
  fwd = fftw_plan_many_dft(3,ngridArray.data(),ndata,
                           ptr,0,istride,ngridN,
                           ptr,0,istride,ngridN,
                           FFTW_FORWARD,FFTW_MEASURE);
  rev = fftw_plan_many_dft(3,ngridArray.data(),ndata,
                           ptr,0,istride,ngridN,
                           ptr,0,istride,ngridN,
                           FFTW_BACKWARD,FFTW_MEASURE);
}

PairIntegrator::~PairIntegrator() {
  fftw_destroy_plan(rev);
  fftw_destroy_plan(fwd);
}

void PairIntegrator::integrate(double q, double scaleTau) {
  double tauSave = tau;
  tau *= scaleTau;
  blitz::Range all(blitz::Range::all());
  // Initialize potential grid.
  for (int i=0; i<ngrid; ++i) {
    double x2 = ((i-ngrid/2)*dr)*((i-ngrid/2)*dr); 
#if NDIM>1
    for (int j=0; j<ngrid; ++j) {
      double y2 = ((j-ngrid/2)*dr)*((j-ngrid/2)*dr); 
#if NDIM>2
      for (int k=0; k<ngrid; ++k) {
        double z2 = ((k-ngrid/2)*dr+q)*((k-ngrid/2)*dr+q);
        vgrid(i,j,k) = pot(sqrt(x2+y2+z2));
      }
#else
      vgrid(i,j) = pot(sqrt(x2+y2));
#endif
    }
#else
    vgrid(i) = pot(sqrt(x2));
#endif
  }
  // Initialize wavefunction.
  psi = 0.0; psi0 = 0.0;
  double a = pow(dr,-NDIM);
  int ndiff = ngrid/2/8/(norder>0?norder:1);
  if (ndiff==0) ndiff=1;
  IVecN1 index1 = ngrid/2;
  IVecN2 index2 = ngrid/2;
  int idata = 0;
  for (int iorder=0; iorder<norder+1; ++iorder) {
    index1(NDIM) = ngrid/2 + iorder*ndiff;
    index2(NDIM+1) = ngrid/2 + iorder*ndiff;
    for (int ioff=0; ioff < ((NDIM>1)?iorder+1:1); ++ioff) {
      index1(0) = idata;
      index2(1) = idata;
#if NDIM>1
      index1(1) = ngrid/2 + ioff*ndiff;
      index2(2) = ngrid/2 + ioff*ndiff;
#endif
      for (int i=0; i<maxiter; ++i) {
        index2(0) = i;
        psi(index2) = a;
      }
      psi0(index1) = a;
      ++idata;
    }
  }
  //Propagate the wavefunction forward in a series of segments.
  for (int iseg=0; iseg<nsegment; ++iseg) {
    std::cout << "Segment " << iseg+1 << " (of " << nsegment
              << ", q=" << q << ", tol=" << tol << "):" << std::endl;
    propagate(tau/nsegment, tol);
  }
  // Calculate free propagator (includes FFT normalization).
  double fftNorm = 1./ngridN;
  fftw_complex *ptr0 = (fftw_complex*)psi0.data();
  Array tFlat(tgrid.data(),IVec1(ngridN),blitz::neverDeleteData);
  CArray2 psi0Flat(psi0.data(),IVec2(ndata,ngridN),blitz::neverDeleteData);
  expFullTtau = exp(-tFlat*tau)*fftNorm;
  fftw_execute_dft(fwd,ptr0,ptr0);
  for (int idata=0; idata<ndata; ++idata) {
    for (int igrid=0; igrid<ngridN; ++igrid) {
      psi0Flat(idata,igrid) *= expFullTtau(igrid);
    }
  }
  fftw_execute_dft(rev,ptr0,ptr0);
  // Calculate arrays of G, G0, s2, and z.
  idata = 0;
  Vec r1=0., r2=0.;
  index2(0) = 0;
  for (int iorder=0; iorder<norder+1; ++iorder) {
    r1(NDIM-1) = q + iorder*ndiff*dr;
    r2(NDIM-1) = q - iorder*ndiff*dr;
    index1(NDIM) = ngrid/2 - iorder*ndiff;
    index2(NDIM+1) = ngrid/2 - iorder*ndiff;
    for (int ioff=0; ioff < ((NDIM>1)?iorder+1:1); ++ioff) {
#if NDIM>1
      r1(0) = ioff*ndiff*dr;
      r2(0) = -ioff*ndiff*dr;
      index1(1) = ngrid/2 - ioff*ndiff;
      index2(2) = ngrid/2 - ioff*ndiff;
#endif
      index1(0) = idata;
      index2(1) = idata;
      g(idata) = psi(index2).real();
      g0(idata) = psi0(index1).real();
      s2(idata) = dot(r1-r2,r1-r2)/(q*q);
      z(idata) = (sqrt(dot(r1,r1)) - sqrt(dot(r2,r2)))/q;
      ++idata;
    }
  }
  //Now find coefficients G(s^2)=sum c_ij s^2j z^2i
  Array2 mat(ndata,ndata);
  IArray ipiv(ndata);
  u = log(g0/g);
  int jdata = 0;
  for (int jorder=0; jorder<norder+1; ++jorder) {
    for (int joff=0; joff < ((NDIM>1)?jorder+1:1); ++joff) {
      for (int i=0; i<ndata; ++i) {
        mat(jdata,i) = pow(s2(i),jorder-joff)*pow(z(i),2*joff);
      }
      ++jdata;
    }
  }
  int info, nrhs=1;
  dgesv_(&ndata,&nrhs,mat.data(),&ndata,ipiv.data(),u.data(),&ndata,&info);
  //Reset tau value.
  tau = tauSave;
}

void PairIntegrator::propagate(double segTau, double tol) {
  blitz::Range all(blitz::Range::all());
  int checkiter = 0;
  for (int iter=0; iter<maxiter; ++iter) {
    blitz::Range range(blitz::Range(0,iter));
    IVec3 shape(maxiter,ndata,ngridN);
    CArray3 psiFlat(psi.data(),shape,blitz::neverDeleteData);
    fftw_complex* ptr = (fftw_complex*)(&psi(iter,0,0));
    Array vFlat(vgrid.data(),IVec1(ngridN),blitz::neverDeleteData);
    Array tFlat(tgrid.data(),IVec1(ngridN),blitz::neverDeleteData);
    // Set delta tau.
    double deltaTau=segTau/nstep(iter);
    // Initialize propagators (put FFT normalization in T propagator).
    double fftNorm = 1./ngridN;
    expVtau = exp(-vFlat*deltaTau); 
    expHalfVtau = exp(-0.5*vFlat*deltaTau);
    expTtau = exp(-tFlat*deltaTau)*fftNorm;
    // Advance psi.
    for (int idata=0; idata<ndata; ++idata) {
      for (int igrid=0; igrid<ngridN; ++igrid) {
        psiFlat(iter,idata,igrid) *= expHalfVtau(igrid);
      }
    }
    for (int istep=0; istep<nstep(iter)-1; ++istep) {
      fftw_execute_dft(fwd,ptr,ptr);
      for (int idata=0; idata<ndata; ++idata) {
        for (int igrid=0; igrid<ngridN; ++igrid) {
          psiFlat(iter,idata,igrid) *= expTtau(igrid);
        }
      }
      fftw_execute_dft(rev,ptr,ptr);
      for (int idata=0; idata<ndata; ++idata) {
        for (int igrid=0; igrid<ngridN; ++igrid) {
          psiFlat(iter,idata,igrid) *= expVtau(igrid);
        }
      }
    }
    fftw_execute_dft(fwd,ptr,ptr);
    for (int idata=0; idata<ndata; ++idata) {
      for (int igrid=0; igrid<ngridN; ++igrid) {
        psiFlat(iter,idata,igrid) *= expTtau(igrid);
      }
    }
    fftw_execute_dft(rev,ptr,ptr);
    for (int idata=0; idata<ndata; ++idata) {
      for (int igrid=0; igrid<ngridN; ++igrid) {
        psiFlat(iter,idata,igrid) *= expHalfVtau(igrid);
      }
    }
/*
std::ostringstream ss; ss << "psi" << iter << ".dat";
std::ofstream file(ss.str().c_str());
for (int i=0; i<ngrid; ++i) {
  for (int j=0; j<ngrid; ++j) {
    file << psi(iter,0,i,ngrid/2,j).real() << std::endl;
  }
  file << std::endl;
}
file.close();
*/
    // Extrapolate.
    if (iter==maxiter-1 || iter > checkiter) {
      vpolyfit(nstepInv(range), psiFlat, best, err, worka, workc, workd);
      double maxerr=0., maxpsi=0.;
      for (int idata=0; idata<ndata; ++idata) {
        for (int igrid=0; igrid<ngridN; ++igrid) {
          double value = best(idata,igrid).real();
          double error = err(idata,igrid).real();
          if (value>maxpsi) maxpsi=value;
          if (error>maxerr) maxerr=error;
        }
      }
      maxerr /= maxpsi;
      std::cout << "  nstep=" << nstep(iter) 
                << ", maxerr=" << maxerr << std::endl;
      if (iter==maxiter-1 || maxerr < tol) {
        for (int i=0; i<maxiter; ++i) {
          for (int idata=0; idata<ndata; ++idata) {
            for (int igrid=0; igrid<ngridN; ++igrid) {
              psiFlat(i,idata,igrid) = best(idata,igrid);
            }
          }
        }
        break;
      } else {
        checkiter += (int)ceil(0.1*log(maxerr/tol));
      }
    }
  }
}

void PairIntegrator::vpolyfit(const Array &x, const CArray3 &y, CArray2 &y0, 
    CArray2 &diff, CArray2 &a, CArray3 &c, CArray3 &d) {
  int n=x.size();
  IVec3 shape(y.shape());
  int n1=shape[1], n2=shape[2];
  blitz::Range all(blitz::Range::all());
  int ntot = 2*n*shape[1]*shape[2];
  int none = 2*shape[1]*shape[2];
  int one = 1;
  //y0 = y(n-1,all,all);
  //c(blitz::Range(0,n-1),all,all) = y;
  //d(blitz::Range(0,n-1),all,all) = y;
  dcopy_(&none,(double*)&y(n-1,0,0),&one,(double*)y0.data(),&one);
  dcopy_(&ntot,(double*)y.data(),&one,(double*)c.data(),&one);
  dcopy_(&ntot,(double*)y.data(),&one,(double*)d.data(),&one);
  for (int j=1; j<n; ++j) {
    for (int i=0; i<n-j; ++i) {
      double denom=1./(x(i)-x(i+j));
      dcopy_(&none,(double*)&c(i+1,0,0),&one,(double*)a.data(),&one);
      dscal_(&none,&denom,(double*)a.data(),&one);
      denom *= -1;
      daxpy_(&none,&denom,(double*)&d(i,0,0),&one,(double*)a.data(),&one);
      dcopy_(&none,(double*)a.data(),&one,(double*)&c(i,0,0),&one);
      dscal_(&none,&x(i),(double*)&c(i,0,0),&one);
      dcopy_(&none,(double*)a.data(),&one,(double*)&d(i,0,0),&one);
      dscal_(&none,&x(i+j),(double*)&d(i,0,0),&one);
/*      for (int i1=0; i1<n1; ++ i1) {
        for (int i2=0; i2<n2; ++ i2) {
          a(i1,i2) = (c(i+1,i1,i2)-d(i,i1,i2))*denom;
          c(i,i1,i2) = x(i)*a(i1,i2);
          d(i,i1,i2) = x(i+j)*a(i1,i2);
        }
      } */
    }
    //y0 += d(n-j-1,all,all);
    double unity=1.;
    daxpy_(&none,&unity,(double*)&d(n-j-1),&one,(double*)y0.data(),&one);
  }
  diff = d(0,all,all);
}
