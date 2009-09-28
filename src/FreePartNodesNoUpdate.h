// $Id$
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
#ifndef __FreePartNodesNoUpdate_h_
#define __FreePartNodesNoUpdate_h_

#include "DoubleAction.h"
#include <blitz/array.h>
#include <blitz/tinymat.h>
#include <vector>
class SimulationInfo;
class Paths;
class SuperCell;
class Species;
#include "PeriodicGaussian.h"

/** Free particle nodes for new fermion algorithm.
This nodal restriction rejects any paths for which
@f[ \rho_0(R(t),R(t+\beta/2),\beta/2) < 0; \qquad (0\le t < \beta/2). @f]
It does this by returning HUGE_NUMBER for the action value from
getActionDifference and getTotalAction methods.

The free particle slater determinant is a slater determinant of
PeriodicGaussians with width @f$\sigma^2=\beta/2m@f$, 
or @f$\alpha=1/(mT)@f$.

There is a nodal action for paths that come close to the node but
do not cross it. As in Ceperley's article, we approximate this
fermion action as
@f[ U = -\log [ 1-\exp(-2m d_i d_{i-1}/ \tau)], @f]
where @f$\tau@f$ is the link timestep, @f$d_i@f$ is the distance from
bead @f$i@f$ to the node, and @f$m@f$ is is the mass of a fermion.
We estimate the distance to the node as the inverse log gradient,
@f$d_i\approx\rho_T/|\nabla_{R_i}\rho_T|@f$.
This gives a nodal contribution to the beta derivative of the action,
@f[\dot{U_i}=\frac{X_{i-1,i} e^{-X_{i-1,i}}}{\tau(1-e^{-X_{i-1,i}})}, @f]
where @f$X_{i-1,i}=2md_{i-1}d_i/\tau@f$. Here we have ignored the
@f$\tau@f$ dependence of the nodal distance, since it is negligible
when @f$\tau\ll\beta@f$.

For the virial estimator, we need the nodal force, defined as the
negative gradient of the nodal action for coordinates at slice i.
From the chain rule we get
@f[ F = -\nabla_{R_i} U = 
\left[ X_{i+1,i}\frac{e^{-X_{i+1,i}}}{1-e^{-X_{i+1,i}}}
     + X_{i,i-1}\frac{e^{-X_{i,i-1}}}{1-e^{-X_{i,i-1}}} \right]
\frac{\nabla d_i}{d_i}. @f]
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

Before sampling, once we know the section, a call to initialize
evalutates and stores inverse slater matricies and the distances
to the nodes for each slice in the section, including the end
slices. The evaluateUpdate method calculates the new slater matrices
and (on the last level) the new distances to the nodes. Evaluate
update returns the difference in nodal action (HUGE_NUMBER if a
node was crossed). If the move is accepted, acceptLastMove updates
the stored matrices and nodal distances.

@version $Revision$
@author John Shumway */
class FreePartNodesNoUpdate : public DoubleAction {
public:
  /// Constants and typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::Array<Mat,2> MMatrix;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  FreePartNodesNoUpdate(const SimulationInfo&, const Species&,
                        const double temperature, const int maxlevel=10);
  /// Destructor.
  ~FreePartNodesNoUpdate();
  /// Calculate the difference in action.
  virtual double getActionDifference(const DoubleMLSampler&,int level);
  virtual double getActionDifference(const DoubleDisplaceMoveSampler&,
                                     const int nMoving) {std :: cout << "????? return 0"<<std ::endl; return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  /// Evaluate the slater determinant for a slice.
  double evaluate(const VArray& r1, const VArray& r2, const int islice) const;
  /// Evaluate the slater determinant for a slice.
  double evaluateDistance(const int islice) const;
  // Evaluate gradient of log of Slater determinant.
  //virtual double evaluateGrad(const VArray& r1, const VArray& r2, const int);
  /// Evaluate the slater determinant for all paths slice.
  double evaluate(const Paths&, const VArray& ref, const int nSlice) const;
  /// Initialize for a sampling section.
  virtual void initialize(const DoubleSectionChooser&);
  /// Accept last move.
  virtual void acceptLastMove();
private:
  /// The time step.
  const double tau;
  /// The mass.
  const double mass;
  /// Number of particles of this type of fermion.
  const int npart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Number of slices.
  int nslice;
  /// Storage for the current slice.
  mutable VArray r1,r2;
  /// The slater determinants.
  mutable std::vector<Matrix*> matrix;
  /// Gradient of the slater determinant.
  mutable VMatrix gradmat;
  /// Determinants.
  mutable Array det, newDet;
  /// Distances to the nodes.
  mutable Array dist, newDist;
  /// Matrix diagonalization arrays.
  mutable IArray ipiv;
  const int lwork;
  mutable Array work;
  /// The SuperCell.
  SuperCell& cell;
  /// A periodic gaussian.
  PeriodicGaussian pg;
  /// Flag for checking if this species is being moved.
  bool notMySpecies;
  /// Storage for forces.
  mutable VArray force;
  /// Storage for first derivatives needed for forces.
  mutable VArray gradArray;
  mutable VMatrix gradMatrix;
  /// Storage for second derivatives needed for forces.
  mutable MMatrix grad2Matrix;
};
#endif
