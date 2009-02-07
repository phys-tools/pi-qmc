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
#ifndef __FixedNodeAction_h_
#define __FixedNodeAction_h_

#include "DoubleAction.h"
#include "NodeModel.h"
#include <blitz/array.h>
#include <blitz/tinymat.h>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include <vector>
class SimulationInfo;
class Paths;
class Species;

/** Action for fixed node approximation.
This nodal restriction rejects any paths for which
@f[ \rho_0(R(t),R(t+\beta/2),\beta/2) < 0; \qquad (0\le t < \beta/2). @f]
The function @f$\rho_0@f$ is a real-valued scalar defined in the
NodeModel. The rejection is implemented 
by returning HUGE_NUMBER for the action value from
getActionDifference and getTotalAction methods.

There is a nodal action for paths that come close to the node but
do not cross it. As in Ceperley's article, we approximate this
fermion action as
@f[ U = -\log [ 1-\exp(-2m d_i d_{i-1}/ \tau)], @f]
where @f$\tau@f$ is the link timestep, @f$d_i@f$ is the distance from
bead @f$i@f$ to the node, and @f$m@f$ is is the mass of a fermion.
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

Before sampling, once we know the section, a call to initialize
evalutates and stores the NodeModel values and the distances
to the nodes for each slice in the section, including the end
slices. The evaluateUpdate method calculates the new NodeModel values
and (on the last level) the new distances to the nodes. Evaluate
update returns the difference in nodal action (HUGE_NUMBER if a
node was crossed). If the move is accepted, acceptLastMove updates
the NodeModel values and nodal distances.

@version $Revision$
@author John Shumway */
class FixedNodeAction : public DoubleAction {
public:
  /// Constants and typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef tvmet::Matrix<double,NDIM,NDIM> Mat;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<double,2> Array2;
  typedef blitz::Array<double,3> Array3;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::Array<Mat,2> MMatrix;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  FixedNodeAction(const SimulationInfo&, const Species&, NodeModel*,
                  const bool withNodalAction=true,
                  const bool useDistDerivative=false, const int maxlevel=12);
  /// Destructor.
  ~FixedNodeAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const DoubleMLSampler&,int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  /// Initialize for a sampling section.
  virtual void initialize(const DoubleSectionChooser&);
  /// Accept last move.
  virtual void acceptLastMove();
private:
  /// The time step.
  const double tau;
  /// Total number of particles.
  const int npart;
  /// Number of particles of this type of fermion.
  const int nSpeciesPart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Number of slices.
  int nslice;
  /// Storage for the current slice.
  mutable VArray r1,r2;
  /// Flag for checking if this species is being moved.
  bool notMySpecies;
  /// Determinants.
  mutable Array dmValue, newDMValue;
  /// Distances to the nodes.
  mutable Array3 dist, newDist;
  /// Array to cache forces for getBeadAction.
  mutable VArray force;
  mutable VMatrix gradd1, gradd2;
  /// Arrays to store distances for getBeadAction.
  mutable Array dim1, dip1, di1, dim2, dip2, di2, 
                dotdim1, dotdi1, dotdim2, dotdi2; 
  /// The NodeModel.
  NodeModel* nodeModel;
  /// The nodeModel's MatrixUpdate object.
  NodeModel::MatrixUpdate* matrixUpdateObj;
  /// Flag for including the nodal action when near the node.
  bool withNodalAction;
  /// Flag for calculating time derivative of the distance to the node.
  bool useDistDerivative;
  int nerror;
};
#endif
