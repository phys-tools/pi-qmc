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
#ifndef __SpinPhase_h_
#define __SpinPhase_h_

#include "SpinPhaseModel.h"
#include <vector>
#include <complex>
#include <blitz/array.h>
class SimulationInfo;
class Species;

/** Class for Spin phase action in the 4-D Spin simulation.
*
* A spin-1/2 particle in a magnetic field has energy
* @f[ H = g\mu_B sB\; \vec{\sigma}\cdot\hat{B} @f]
* where the magnetic field is @f$\vec{B}=B\hat{B}@f$.
*
* We represent the up and down spin states with complex functions,
* @f[ \chi_\uparrow(\vec{s}) = \pi^{-1}(s_0+is_3)e^{-s^2/2}
*   \qquad  \chi_\downarrow(\vec{s}) = \pi^{-1}(-s_2-is_1)e^{-s^2/2}.
* @f]
* The thermal density matrix is 
* @f$ \hat{\rho} = \frac{1}{2}(\hat{\mathbb 1} + \hat{B}\cdot\vec{\sigma}
*     \tanh x), @f$ where @f$x=\beta g \mu_B s B@f$.
*
* @f[ \rho(\vec {s}, \vec {s'}) = f(\vec {s} , \vec {s'})[\vec {s}\cdot \vec{s'} 
*   +i(s_{0}s'_{3}-s'_{0}s_{3} + s_{1}s'_{2}-s'_{1}s_{2})]   @f]
*  where @f[ f(\vec {s} , \vec {s'}) = \frac{1}{16\pi^{2}} e^{-\frac{(s^{2}+s'^{2})}{2}} @f] 
*
* The derivative of the density matrix are 
* @f[ \nabla \rho(\vec{s},\vec{s'})=\vec{s'}+i(-s'_{3}\hat{0}-s'_{2}\hat{1}
*              +s'_{1}\hat{2}  +s'_{0}\hat{3})  @f] 
* @f[ \nabla' \rho(\vec{s},\vec{s'})=\vec{s}+i(s_{3}\hat{0'}+s_{2}\hat{1'}
*     		+s_{1}\hat{2'}  +s_{0}\hat{3'})  @f] 
*  where @f$ \hat{0}, \hat{1}, \hat{2},  \hat{3} @f$ are the unit vector of @f$ \vec{s} @f$.
* 
** @author Daejin Shin and John Shumway  */
class SpinPhase : public SpinPhaseModel {
public:
  typedef std::complex<double> Complex;
  typedef blitz::TinyVector<double,4> SVec;
  typedef blitz::TinyVector<Complex,4> C4Vec;
  typedef blitz::TinyVector<Complex,NDIM> CVec;
  typedef blitz::Array<Complex,2> Matrix;
  typedef blitz::Array<CVec,2> VMatrix;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Complex,1> CArray;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  SpinPhase(const SimulationInfo&, const Species&, const double t,
	    const double bx, const double by, const double bz,
	    const double gmubs, const int maxlevel);
  /// Virtual destructor.
  virtual ~SpinPhase();
  /// Evaluate the phase, gradient of of the phase, and vector potential.
  virtual void evaluate(const VArray&, const VArray&,
                        const SArray&, const SArray&, const int islice);
  /// Evaluate a one-particle density matrix.
  void evalOrbit(const Vec& r1, const Vec& r2, const SVec& s1, const SVec& s2,
                 Complex &rho, Vec &gradrho1, Vec &gradrho2,
                 C4Vec &sgradrho1, C4Vec &sgradrho2);
private:
  /// The time step.
  const double tau;
  /// The charge.
  const double charge;
  /// The temperature.
  const double temperature;
  /// The magnetic field strength.
  const double bx, by, bz;
  /// The constant @f$ g\mu_B s @f$.
  const double gmubs;
  /// The value of @f$ \tanh \beta g\mu_B B s @f$.
  const double tanhx, tanhy, tanhz;
  /// Number of particles of this type of fermion.
  const int npart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Number of slices.
  int nslice;
  /// Inverse slater matrices, stored at each slice.
  std::vector<Matrix*> matrix; 
  /// Temporary storage for gadient of the slater matrix.
  VMatrix gradmat1, gradmat2; 
  /// Matrix diagonalization arrays.
  IArray ipiv;
  const int lwork;
  CArray work;
  /// Hyperbolic constants.
  const double f;
  /// Pi.
  static const double pi;
};
#endif
