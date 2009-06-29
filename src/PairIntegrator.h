// $Id: PairIntegrator.h 50 2009-05-11 20:05:24Z john.shumwayjr $
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
#ifndef __PairIntegrator_h_
#define __PairIntegrator_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
class Species;
class PairPotential;
#include <blitz/array.h>
#include <fftw3.h>

/// Class for integrating pair potential for finite timestep.
/// @author John Shumway. 
class PairIntegrator {
public:
  typedef std::complex<double> Complex;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<double,2> Array2;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Complex,1> CArray;
  typedef blitz::Array<Complex,2> CArray2;
  typedef blitz::Array<Complex,3> CArray3;
  typedef blitz::Array<double,NDIM> ArrayN;
  typedef blitz::Array<Complex,NDIM> CArrayN;
  typedef blitz::Array<Complex,NDIM+1> CArrayN1;
  typedef blitz::Array<Complex,NDIM+2> CArrayN2;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,1> IVec1;
  typedef blitz::TinyVector<int,2> IVec2;
  typedef blitz::TinyVector<int,3> IVec3;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::TinyVector<int,NDIM+1> IVecN1;
  typedef blitz::TinyVector<int,NDIM+2> IVecN2;
  PairIntegrator(double tau, double mu, double dr, int norder, int maxiter,
                 const PairPotential &pot);
  ~PairIntegrator();
  void integrate(double q);
  void propagate(double segTau, double tol);
  void vpolyfit(const Array &x, const CArray3 &y, CArray2 &y0, 
                CArray2 &diff, CArray2 &a, CArray3 &c, CArray3 &d);
  const Array& getU() {return u;}
private:
  const double tau;
  const double mu;
  double dr;
  const int norder;
  const int maxiter;
  const int ndata;
  int ngrid, ngridN;
  IArray nstep;
  Array nstepInv;
  ArrayN vgrid;
  ArrayN tgrid;
  Array expVtau, expHalfVtau, expTtau, expFullTtau;
  CArrayN2 psi;
  CArrayN1 psi0;
  CArray3 workc, workd;
  CArray2 best, err, worka;
  Array g,g0,s2,z,u;
  fftw_plan fwd, rev;
  static const double PI=3.1415926535897931;
  const PairPotential &pot;
};
#endif
