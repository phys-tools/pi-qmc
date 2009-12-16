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
#ifndef __CoulombAction_h_
#define __CoulombAction_h_
class MultiLevelSampler;
class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
class EwaldSum;
#include "Action.h"
#include "PairAction.h"
#include <vector>

/** Class for calculating coulomb action between particles.
  * @version $Revision$
  * @bug getEField method is not very accurate or efficient.
  * @author John Shumway. */
class CoulombAction : public Action, public PairAction::EmpiricalPairAction {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<Vec,1> VArray1;
  /// Constructor by providing the timestep tau.
  CoulombAction(const double epsilon, const SimulationInfo&, const int norder,
                double rmin, double rmax, int ngpts, const bool dumpFiles,
                bool useEwald, int ewaldNDim, double ewaldRcut, 
                double ewaldKcut, double screenDist, const double kappa, const int nimages);
  /// Virtual destructor.
  virtual ~CoulombAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the difference in action.
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int nslice);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  /// Calculate the efield at a bead.
  double getEField(const Paths&, const int ipart, const int islice) const;
  /// Calculate the total action for a slice.
  double getAction(const Paths&, const int islice) const;
  /// EmpiricalPairAction method.
  virtual double u(double r, int order) const;
  /// EmpiricalPairAction method.
  virtual double utau(double r, int order) const;
private:
   /// The dielectric constant.
  const double epsilon;
  /// The timestep.
  mutable double tau;
  /// Number of particles.
  int npart;
  /// The gridded pair-actions.
  std::vector<PairAction*> pairActionArray;
  /// Charge for EmpiricalPairAction method.
  double q1q2;
  /// Reduced mass for EmpiricalPairAction method.
  double mu;
  /// Displacement variable for EmpricalPairAction method.
  double displace2;
  /// Screening distance (zero for no screening).
  double screenDist;
  /// Ewald sums.
  EwaldSum *ewaldSum;
  /// Buffer for positions in long range ewald sum.
  mutable VArray1 rewald;
};
#endif
