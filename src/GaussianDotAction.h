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
#ifndef __GaussianDotAction_h_
#define __GaussianDotAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
#include "Action.h"
#include <blitz/array.h>

/** Class for calculating an attractive Gaussian potenial.
  * @version $Revision$
  * @author John Shumway. */
class GaussianDotAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor.
  GaussianDotAction(const double v0, const double alpha, const Vec &center,
                    const SimulationInfo&);
  /// Virtual destructor.
  virtual ~GaussianDotAction() {}
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
  /// The center.
  Vec center;
};
#endif
