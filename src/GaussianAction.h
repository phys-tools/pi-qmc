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
#ifndef __GaussianAction_h_
#define __GaussianAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
#include "Action.h"
#include <blitz/array.h>

/** Class for calculating coulomb action between particles.
  * @version $Revision$
  * Gives the same interaction between all pairs of particles.
  * @bug Assumes gaussians are small enough to neglect periodic images.
  * @author John Shumway. */
class GaussianAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<int,1> IArray;
  /// Constructor by providing the timestep tau.
  GaussianAction(const double v0, const double alpha, const SimulationInfo&);
  /// Virtual destructor.
  virtual ~GaussianAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
 virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The strength of the interactions.
  const double v0;
  /// The coefficient on the gaussian.
  const double alpha;
  /// The timestep.
  const double tau;
};
#endif
