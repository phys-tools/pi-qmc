//$Id$
/*  Copyright (C) 2004-2008 John B. Shumway, Jr.

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
#ifndef __WireNodes_h_
#define __WireNodes_h_

#include "NodeModel.h"
#include <vector>
#include <cstdlib>
#include <blitz/tinymat.h>
class PeriodicGaussian;
class SimulationInfo;
class Species;
class SuperCell;

/** Free particle nodes.
We define the node model as a slater determinant of density matricies
of single, free particles,
@f[\rho_T(R,R')=\operatorname{det}|\rho(r_j,r_i)|,@f]
where 
@f[\rho(r_j,r_i)=\exp\left(-\frac{m|r_j-r_i|^2}{2\tau}\right).@f]
For periodic boundary conditions, we use must use a PeriodicGaussian.
Note that the normalization factor is not needed for a NodeModel.

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
where the indicies refer to particles.
@author John Shumway */
class WireNodes : public NodeModel {
public:
  typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::Array<Mat,2> MMatrix;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<int,2> IArray2;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Class for matrix updates.
  class MatrixUpdate : public NodeModel::MatrixUpdate {
  public:
    MatrixUpdate(int maxMovers, int maxlevel, int npart,
                 std::vector<Matrix*>& matrix, const WireNodes &wireNodes);
    virtual double evaluateChange(const DoubleMLSampler&, int islice);
    virtual void evaluateNewInverse(const int islice);
    virtual void evaluateNewDistance(const VArray &r1, const VArray &r2,
         const int islice, Array &d1, Array &d2);
    virtual void acceptLastMove(int nslice);
  private:
    const WireNodes& wireNodes;
    const int maxMovers;
    const int npart;
    /// Temporary storage for inverse slater matricies.
    mutable std::vector<Matrix*> newMatrix;
    /// Temporary storage for new slater matrix columns.
    mutable std::vector<Matrix*> phi;
    /// Temporary storage for new slater matrix products.
    mutable std::vector<Matrix*> bvec;
    /// Temporary storage for small (moving) slater determinant.
    mutable Matrix smallDet;
    /// Reference to the inverse slater matrices.
    std::vector<Matrix*>& matrix;
    /// Matrix diagonalization arrays.
    mutable IArray ipiv;
    const int lwork;
    mutable Array work;
    /// Flag for updating inverse.
    bool isNewMatrixUpdated;
    /// Reference to moving particles.
    IArray index1,mindex1;
    int nMoving;
  };
  /// Constructor.
  WireNodes(const SimulationInfo&, const Species&, const double omega,
            const double temperature, const Vec& center, const int maxlevel,
            const bool useUpdates, const int maxMovers);
  /// Virtual destructor.
  virtual ~WireNodes();
  /// Evaluate the density matrix function, returning the value.
  virtual DetWithFlag evaluate(const VArray &r1, const VArray &r2, 
                               const int islice, bool scaleMagnitude);
  /// Evaluate distance to the node in units of @f$ \sqrt{\tau/2m}@f$.
  /// Assumes that evaluate has already been called on the slice.
  virtual void evaluateDistance(const VArray &r1, const VArray &r2,
                                const int islice, Array &d1, Array &d2);
  /// Evaluate the time-derivative of the distance to the 
  /// node in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateDotDistance(const VArray &r1, const VArray &r2,
                                   const int islice, Array &d1, Array &d2);
  /// Evaluate gradient of log of the distance to the node
  /// in units of @f$ \sqrt{\tau/2m}@f$.
  virtual void evaluateGradLogDist(const VArray &r1, const VArray &r2,
           const int islice, VMatrix &gradd1, VMatrix &gradd2,
                             const Array &d1, const Array &d2);
private:
  /// The time step.
  double tau;
  /// The mass.
  const double mass;
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
  /// The SuperCell.
  SuperCell& cell;
  /// A periodic gaussian.
  PeriodicGaussian *pg, *pgp, *pgm;
  /// The center of the harmonic potential.
  const Vec center;
  /// The wire confinement strength.
  const double omega; 
  /// More parameters.
  const double coshwt, sinhwt, c;
  /// Flag for checking if this species is being moved.
  bool notMySpecies;
  /// Storage for first derivatives needed for forces.
  mutable VArray gradArray1, gradArray2;
  //mutable VArray gradArray;
  //mutable VMatrix gradMatrix, mat2;
  /// Storage for second derivatives needed for forces.
  //mutable MMatrix grad2Matrix;
  /// Step for finite difference calculation of time derivative.
  static const double EPSILON;
  /// Storage for calculating derivatives.
  Array temp1, temp2;
  /// Storage for calculating dominant contribution to determinant.
  Matrix uarray;
  IArray2 kindex;
  IArray kwork;
  int nerror;
  double scale;
};
#endif
