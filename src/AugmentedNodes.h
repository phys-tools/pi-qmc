//$Id: AugmentedNodes.h 52 2009-05-13 18:51:51Z john.shumwayjr $
/*  Copyright (C) 2009,2010 John B. Shumway, Jr.

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
#ifndef __AugmentedNodes_h_
#define __AugmentedNodes_h_

#include "NodeModel.h"
#include <vector>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
class PeriodicGaussian;
class SimulationInfo;
class Species;
class SuperCell;

/** Augmented free particle nodes.
We define the node model as a slater determinant of density matricies
of single particles,
@f[\rho_T(R,R')=\operatorname{det}|\rho(r_j,r_i)|.@f]
The single particle density matrices are taken to be a mixture
of free particle density matricies plus orbital density matrices,
@f[\rho(r_j,r_i)=n_0 \exp\left(-\frac{m|r_j-r_i|^2}{2\tau}\right)
+ w_i\sum_k \psi(r_j)\psi^*(r_i),@f]
where @f$n_0@f$ is the density of the free particle
density matrix and the @f$w_k@f$ are the weights of the atomic orbitals.
For periodic boundary conditions, we use must use a PeriodicGaussian in 
the free particle density matrix; we assume the atomic orbitals are
smaller than the periodic box.

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

If this works we will want to rewrite for efficiency!!!
@author John Shumway */
class AugmentedNodes : public NodeModel {
public:
  typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<int,2> IArray2;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Class for matrix updates.
  class MatrixUpdate : public NodeModel::MatrixUpdate {
  public:
    MatrixUpdate(int maxMovers, int maxlevel, int npart,
       std::vector<Matrix*>& matrix, const AugmentedNodes &fpnodes);
    virtual double evaluateChange(const DoubleMLSampler&, int islice);
    virtual void evaluateNewInverse(const int islice);
    virtual void evaluateNewDistance(const VArray &r1, const VArray &r2,
      const int islice, Array &d1, Array &d2);
    virtual void acceptLastMove(int nslice);
  private:
    const AugmentedNodes& fpNodes;
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
  /// Classes for atomic orbital density matricies.
  class AtomicOrbitalDM {
  public:
    AtomicOrbitalDM(int nuclearIndex, double weight);
    const double weight;
    const int nuclearIndex;
    struct ValueAndGradient {
      double value, gradr1, gradr2, gradcostheta;
    };
    virtual ValueAndGradient 
    operator()(double r1, double r2, double costheta) const=0;
  };
  class Atomic1sDM : public AtomicOrbitalDM {
  public:
    Atomic1sDM(double Z, int nuclearIndex, double weight);
    virtual ValueAndGradient 
    operator()(double r1, double r2, double costheta) const;
    const double Z;
  };
  class Atomic2sDM : public AtomicOrbitalDM {
  public:
    Atomic2sDM(double Z, int nuclearIndex, double weight);
    virtual ValueAndGradient 
    operator()(double r1, double r2, double costheta) const;
    const double Z;
  };
  class Atomic2pDM : public AtomicOrbitalDM {
  public:
    Atomic2pDM(double Z, int nuclearIndex, double weight);
    virtual ValueAndGradient 
    operator()(double r1, double r2, double costheta) const;
    const double Z;
  };
  /// Constructor.
  AugmentedNodes(const SimulationInfo&, const Species&, const Species&, 
    const double temperature, const int maxlevel, 
    const bool useUpdates, const int maxMovers, double density, 
    const std::vector<const AtomicOrbitalDM*>&, bool useHungarian);
  /// Virtual destructor.
  virtual ~AugmentedNodes();
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
  /// Returns true if action depends on other particle coordinates.
  virtual bool dependsOnOtherParticles() {return true;}
private:
  /// The time step.
  double tau;
  /// The mass.
  const double mass;
  const double mass2;
  /// Number of particles of this type of fermion.
  const int npart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  const int npart2;
  const int kfirst;
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
  std::vector<PeriodicGaussian*> pg, pgp, pgm;
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
  /// Flag for number of bad returns from LAPACK calls.
  int nerror;
  IArray kindex2;
  Matrix kmat;
  /// Flag to use Hungarian method to remove largest contribution to det.
  const bool useHungarian;
  /// Scale factor to avoid overflow or underflow.
  double scale;
  /// Constant.
  static const double PI;
  /// Density of free particle density matrix.
  const double density;
  const std::vector<const AtomicOrbitalDM*> orbitals;
};
#endif
