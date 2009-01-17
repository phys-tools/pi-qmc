//$Id: SHORealNodes.h,v 1.5 2007/10/26 07:04:26 jshumwa Exp $
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
#ifndef __SHORealNodes_h_
#define __SHORealNodes_h_

#include "NodeModel.h"
#include <vector>
#include <complex>
#include <blitz/array.h>
class SimulationInfo;
class Species;

/** Simple hamonic oscillator Real part of the density matrix. 
 ** We can write the density matrix for a harmonic oscillator 
    in a constant magnetic field along with z direction. 
    The Hamiltonian is given by @f[
    H=\frac{1}{m}(\vec{p}+\frac{e}{c}\vec{A})^2 +\frac{1}{2}mw_0^2\vec r^{2} @f]
    The vector potential is @f$ A=(-\frac{1}{2}B_0 y,\frac{1}{2}B_0 x,0 )@f$
    From K. Yonei [ J. Phys.A:Math.Gen.22,(1989) 2415 ], 
    the density matrix can be written as 
    @f[<x|exp(-\beta H)|x'><y|exp(-\beta H)|y'> =  
    \dfrac{mw'}{2\pi sinh\beta w'}exp\Big(-a_1[(x+x')^2 +(y+y')^2]-a_2[(x-x')^2 +(y-y')^2]
              -ib(xy'-x'y) \Big) @f]
    and for z direction,
    @f[ <z|exp(-\beta H)|z'>= \sqrt{\dfrac{mw_0}{2\pi sinh(\tau w')}}
    exp \Big(-[c_1 (z+z')^2 +c_2 (z-z')^2]\Big). @f] 
    ,where @f$ w'=\sqrt{w_0^{2}+w_{c}^{2}} @f$, and  @f$w_c=\dfrac{eB_0}{2m}@f$ 
    The coefficients are the follows.
   @f[  a_1=\dfrac{mw'}{4sinh(\tau w')}[cosh(\tau w') -cosh(\tau w_c)]  @f]
   @f[  a_2=\dfrac{mw'}{4sinh(\tau w')}[cosh(\tau w') +cosh(\tau w_c)]  @f]
   @f[  b=\dfrac{mw'sinh(\tau w_c)}{sinh(\tau w')}                   @f]
   @f[  c_1= \frac{mw_0}{4}tanh(\frac{\tau w_0}{2})                 @f]
   @f[  c_2= \frac{mw_0}{4}coth(\frac{\tau w_0}{2})                 @f]
   We need the x derivative of the thermal density matrix without the imaginary term.
   @f[ \nabla_{x}\rho(x,x';\tau)= [-2a_1((x+x')+(y+y')) -2a_2((x-x')+(y-y'))] 
         exp\Big(-a_1[(x+x')^2 +(y+y')^2]-a_2[(x-x')^2 +(y-y')^2] \Big). @f]
    Suppose we have the density matrix of complex funtion, 
    we can write @f$ \rho = |\rho|e^{i\phi}.  @f$ 
    We need the derivative of density matrix.
 @f[\nabla\phi=Im \left(\dfrac{\nabla\rho}{\rho} \right) @f]
* @author Daejin Shin */
class SHORealNodes : public NodeModel {
public:
  typedef std::complex<double> Complex;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<Complex,NDIM> CVec;
  typedef blitz::Array<Complex,2> Matrix;
  typedef blitz::Array<CVec,2> VMatrix;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<Complex,1> CArray;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  SHORealNodes(const SimulationInfo&, const Species&, const double omega,
           const double temperature, const double b, const int maxlevel);
  /// Virtual destructor.
  virtual ~SHORealNodes();
   /// Evaluate the density matrix function, returning the value.
  virtual double evaluate(const VArray &r1, const VArray &r2, 
                          const int islice);
  /// Evaluate distance to the node in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateDistance(const VArray &r1, const VArray &r2,
                                const int islice, Array &d1, Array &d2);
  /// Evaluate gradient of log of the distance to the node
  /// in units of @f$ \sqrt{\tau/2m}@f$.
  //virtual void evaluateGradLogDist(const VArray &r1, const VArray &r2,
  //                 const int islice, VArray &force, const double distance);
  virtual void evaluateGradLogDist(const VArray &r1, const VArray &r2,
           const int islice, VMatrix &gradd1, VMatrix &gradd2,
                             const Array &d1, const Array &d2)=0;


private:
  /// The time step.
  const double tau;
  /// The temperature.
  const double temperature;
  /// The mass.
  const double mass;
  /// The charge.
  const double charge;
  /// The harmonic confinement frequency.
  const double omega;
  /// The magnetic field.
  const double b;
  /// Number of particles of this type of fermion.
  const int npart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Number of slices.
  int nslice;
  /// Inverse slater matrices, stored at each slice.
  std::vector<Matrix*> matrix; 
  /// Temporary storage for gadient of the slater matrix.
  VMatrix gradmat1, gradmat2, gradmattau; 
  /// Matrix diagonalization arrays.
  mutable IArray ipiv;
  const int lwork;
  mutable CArray work;
  CArray det;
  /// The cyclotron frequency.
  const double omegac;
  /// The effective frequency.
  const double omega1;
  /// Hyperbolic constants.
  const double sinh1,cosh1,sinhc,coshc;
  /// More constants.
  const double a1,a2,b1,a1tau,a2tau,b1tau;
  /// The speed of light.
  static const double c;
  /// Pi.
  static const double pi;
};
#endif
