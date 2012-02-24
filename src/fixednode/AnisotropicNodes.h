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
#ifndef __AnisotropicNodes_h_
#define __AnisotropicNodes_h_

#include "DoubleAction.h"
#include <cstdlib>
#include <blitz/array.h>
#include <vector>
class SimulationInfo;
class Paths;
class SuperCell;
class Species;
class DoubleSectionChooser;
#include "PeriodicGaussian.h"

/** Free particle nodes for new fermion algorithm.
  * @todo Calculate nodal action.
  * @todo Add documentation.
  * @bug Hacked to catch case when this species is not moved.
  * @bug Does note work when  sampling both fermion sections.
  * @version $Revision$
  * @author John Shumway */
class AnisotropicNodes : public DoubleAction {
public:
  /// Constants and typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<double,2> Matrix;
  typedef blitz::Array<Vec,2> VMatrix;
  typedef blitz::ColumnMajorArray<2> ColMajor;
  /// Constructor.
  AnisotropicNodes(const SimulationInfo&, const Species&, const double t,
                    const int maxlevel=9, const int maxmovers=3);
  /// Destructor.
  ~AnisotropicNodes();
  /// Calculate the difference in action.
  virtual double getActionDifference(const SectionSamplerInterface&,int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
          double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  /// Evaluate the slater determinant for a slice.
  double evaluate(const VArray& r1, const VArray& r2, Matrix&);
  /// Find the distances to the nodes.
  void evaluateDistances(const int islice);
  /// Evaluate gradient of log of Slater determinant.
  virtual double evaluateGrad(const int islice);
  /// Evaluate the slater determinant using updates.
  double evaluateUpdate(const VArray& r1, const VArray& r2, const int islice);
  /// Initialize for a sampling section.
  virtual void initialize(const DoubleSectionChooser&);
  /// Accept last move.
  virtual void acceptLastMove();
private:
  /// The temperature.
  const double temperature;
  /// The time step.
  const double tau;
  /// The mass.
  const Vec mass;
  /// The square root of the anisotropic mass.
  const Vec sqrtmass;
  /// Number of particles of this type of fermion.
  const int npart;
  /// Index of first particle of this type of fermion.
  const int ifirst;
  /// Storage for the current slice.
  VArray r1,r2;
  /// The slater matricies.
  std::vector<Matrix*> matrix, newMatrix, newColumn;
  /// The gradients of the new Slater matrix.
  std::vector<VMatrix*> newColGrad;
  /// The distances to the nodes.
  Matrix dist,newDist;
  /// The old and new slater determinants.
  Array det,newDet;
  /// The moving index.
  const IArray* movingIndex;
  /// Matrix diagonalization arrays.
  IArray ipiv;
  const int lwork;
  Array work;
  /// The SuperCell.
  SuperCell& cell;
  /// A periodic gaussian.
  PeriodicGaussian* pg[NDIM];
  /// Flag to prevent updates if this species not moved.
  bool notMySpecies;
const DoubleSectionChooser *sectionchooser;
};
#endif
