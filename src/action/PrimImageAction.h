// $Id: PrimImageAction.h 509 2012-04-04 12:40:28Z john.shumwayjr $
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
#ifndef __PrimImageAction_h_
#define __PrimImageAction_h_
class SectionSamplerInterface;class DisplaceMoveSampler;
class Species;
template <int TDIM> class Beads;
#include "Action.h"
#include <cstdlib>
#include <blitz/array.h>
class SimulationInfo;

/** Class for calculating the primitive action for an image charge at a plane 
 * @author John Shumway & Ian Galbraith. 
 */

class PrimImageAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Construct by providing the timestep tau.
  PrimImageAction(const double d , const double del, const double epsilon, const double epsilonrel, const double vbarr, const SimulationInfo &simInfo,
					const int ndim, const Species&);
  /// Virtual destructor.
  virtual ~PrimImageAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const SectionSamplerInterface&,
                                     int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
private:
  /// The timestep.
  const double tau;
  /// The interaction parameters.
  const double d, del, epsilon, epsilonrel ,vbarr;
  /// The number of dimensions (can be less than NDIM, i.e., to make a wire).
  const int ndim;
  /// The first particle in this interaction.
  const int ifirst;
  /// The number of particles with this interaction.
  const int npart;
	
};
#endif
