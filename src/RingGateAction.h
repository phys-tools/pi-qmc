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
#ifndef __RingGateAction_h_
#define __RingGateAction_h_
class MultiLevelSampler;
class Species;
template <int TDIM> class Beads;
#include "Action.h"
#include <blitz/array.h>
class SimulationInfo;

/* Class for calculating the gate action for a ring
 */

class RingGateAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor
  RingGateAction(const SimulationInfo &simInfo, const double GVolt, const double s, const double theta0, const Species&);
  /// Virtual destructor.
  virtual ~RingGateAction() {}
  /// Calcuate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&, const int level);
  /// Calcuate the difference in action (NOT IMPLEMENTED YET).
  virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice, double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
private:
  /// The timestep.
  const double tau;
  /// The interatction parameters.
  const double GVolt, s, theta0;
  /// The first particle in this interaction.
  const int ifirst;
  /// The number of particles with this interaction.
  const int npart;
}; 
#endif
