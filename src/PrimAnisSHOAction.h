// $Id: PrimAnisSHOAction.h 12 2009-02-07 23:32:51Z john.shumwayjr $
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
#ifndef __PrimAnisSHOAction_h_
#define __PrimAnisSHOAction_h_
class MultiLevelSampler;
class Species;
template <int TDIM> class Beads;
#include "Action.h"
#include <blitz/array.h>
class SimulationInfo;

/** Class for calculating the primitive action for a anisotropic simple harmonic
  * oscillator.
  * @author John Shumway. 
 */

class PrimAnisSHOAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Construct by providing the timestep tau.
  PrimAnisSHOAction(const double a, const double b, const double c ,const SimulationInfo &simInfo,
					const int ndim, const Species&);
  /// Virtual destructor.
  virtual ~PrimAnisSHOAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
private:
  /// The timestep.
  const double tau;
  /// The interaction parameters.
  const double a, b, c;
  /// The number of dimensions (can be less than NDIM, i.e., to make a wire).
  const int ndim;
  /// The first particle in this interaction.
  const int ifirst;
  /// The number of particles with this interaction.
  const int npart;
	
};
#endif
