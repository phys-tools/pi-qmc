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
#ifndef __GateAction_h_
#define __GateAction_h_
class MultiLevelSampler;
class Species;
template <int TDIM> class Beads;
#include "Action.h"
#include <blitz/array.h>
class SimulationInfo;

/*  Class for calculating the square gate action
 */

class GateAction : public Action{
public: 
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor
  GateAction(const SimulationInfo &simInfo, const double GVolt, const double sx,
             const double sy, const double xwidth, const double ywidth,
             const double xoffset, const double yoffset, const Species&);
  /// Virtual destructor.
  virtual ~GateAction() {};
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&, const int level);
  /// calculate the difference in action (NOT IMPLEMENTED YET).
  virtual double getActionDifference(const DisplaceMoveSampler&,
                                     const int nMoving) {return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice, double& u,
                       double& utau, double& ulambda, Vec &fm, Vec &fp) const;
private:
  /// The timestep;
  const double tau;
  /// The interaction parameters.
  const double GVolt, sx, sy, xwidth, ywidth, xoffset, yoffset, normalConst;
  /// The first particle in this interaction.
  const int ifirst;
  /// The number of particles with this interaction.
  const int npart;
};
#endif
