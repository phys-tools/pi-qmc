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
#ifndef __SHONodes_h_
#define __SHONodes_h_

#include "NodeModel.h"
#include <vector>
class SimulationInfo;
class Species;

/** Fermionic nodes for a simple hamonic oscillator.

We define the node model as a slater determinant of thermal 
density matricies of single particles in a harmonic oscillator,
@f[\rho_T(R,R')=\operatorname{det}|\rho(r_j,r_i;\beta)|,@f]
where 
@f[\rho(r_j,r_i)
=\exp\left(-\frac{m|r_j-r_i|^2}{2\beta}\right).@f]

We estimate the distance to the node as the inverse log gradient,
@f$d_i\approx\rho_T/|\nabla_{R_i}\rho_T|@f$.

The gradient of the distance is given by first and second derivatives
of the trial density matrix,
@f[
\frac{\nabla_i d}{d}=
\frac{\nabla_i \rho}{\rho}
-d^2 \sum_j
\frac{\overleftarrow{\vec{\nabla}}_{i,j}\rho}{\rho}\cdot
\frac{\vec{\nabla}_j\rho}{\rho}
-d^2 \sum_j
\frac{\overleftarrow{\vec{\nabla}}'_{i,j}\rho}{\rho}\cdot
\frac{\vec{\nabla}'_j\rho}{\rho}.
@f]

  
@f[ U(x,x';\tau) = \frac{1}{2}\log\frac{2\pi\sinh\omega\tau}{m\omega}
  +\frac{m\omega[(x^2+{x'}^2)\cosh\omega\tau-2xx']}{2\sinh\omega\tau}
  -\frac{1}{2}\log 2\pi\tau/m -\frac{m(x-x')^2}{2\tau}. @f]
We also need the time derivative,
@f[ \dot{U}(x,x';\tau) =
   \frac{\omega}{2}\coth\omega\tau
  +\frac{m\omega^2}{2}\left[(x^2+{x'}^2)(1-\coth^2\omega\tau)
  +\frac{2xx'\cosh\omega\tau}{\sinh^2\omega\tau}\right]
  +\frac{m(x-x')^2}{2\tau^2}-\frac{1}{2\tau}. @f]

For the virial estimator, we neee to calculate forces.
@f[ \nabla_{x}U(x,x';\tau)=
   \frac{m\omega[x\cosh\omega\tau-x']}{\sinh\omega\tau}
   -\frac{m(x-x')}{\tau} @f] 
@author John Shumway and Daejin Shin. */
class SHONodes : public NodeModel {
public:
  typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::Array<Mat,2> MMatrix;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  SHONodes(const SimulationInfo&, const Species&, const double omega,
                    const double temperature, const int maxlevel);
  /// Virtual destructor.
  virtual ~SHONodes();
  /// Evaluate the density matrix function, returning the value.
  virtual DetWithFlag evaluate(const VArray &r1, const VArray &r2, 
                               const int islice, bool scaleMagnitude);
  /// Evaluate distance to the node in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateDistance(const VArray &r1, const VArray &r2,
                                const int islice, Array &d1, Array &d2);
  /// Evaluate gradient of log of the distance to the node
  /// in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateGradLogDist(const VArray &r1, const VArray &r2,
           const int islice, VMatrix &gradd1, VMatrix &gradd2,
                             const Array &d1, const Array &d2);
private:
  /// The time step.
  const double tau;
  /// The temperature.
  const double temperature;
  /// The mass.
  const double mass;
  /// The harmonic confinement frequency.
  const double omega;
  /// More parameters.
  const double coshwt, sinhwt, c;
  /// Number of particles of this type of fermion.
  const int npart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Number of slices.
  int nslice;
  /// The inverse slater matricies.
  mutable std::vector<Matrix*> matrix;
  /// Gradient of the slater determinant.
  mutable VMatrix gradmat;
  /// Matrix diagonalization arrays.
  mutable IArray ipiv;
  const int lwork;
  mutable Array work;
  /// Flag for checking if this species is being moved.
  bool notMySpecies;
  /// Storage for first derivatives needed for forces.
  mutable VArray gradArray;
  mutable VMatrix gradMatrix, mat2;
  /// Storage for second derivatives needed for forces.
  mutable MMatrix grad2Matrix;
  int nerror;
};
#endif
