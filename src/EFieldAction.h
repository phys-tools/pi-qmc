// $Id:
/*  Copyright (C) 2004-2006,2010 John B. Shumway, Jr.

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
#ifndef __EFieldAction_h
#define __EFieldAction_h
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
class SuperCell;
#include "Action.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for getting action from stepwise electric field.
  * 
  * @bug Hard coded for NDIM=3.
  * @version $Revision:
  * @author John Shumway and Matthew Harowitz. */
class EFieldAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<double, 1> Array;
  typedef blitz::Array<int, 1> IArray;
  typedef blitz::TinyVector<double, NDIM> Vec;
  /// Construct by providing SimulationInfo, field strength in
  /// atomic units @f$\mathrm{Ha}/a_0 e@f$, physical region center and width,
  /// and direction index.
  EFieldAction(const SimulationInfo& simInfo, double strength,
     double center, double width, int idir);
  /// Virtual destructor.
  virtual ~EFieldAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
   virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
/// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
          double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
  /// The charges of the particles.
  Array q;
  /// The strength of the electric field.
  const double strength;
  /// The center of the physical region.
  Vec center;
  /// The half width of the physical region.
  const double halfwidth;
  /// The component of the electric field (0->x, 1->y, 2->z)
  const int idir;
  /// The supercell
  SuperCell& cell;
  /// The electric potential
  inline double v(const double z) const;
  /// Parameters for v method.
  const double slopeIn, slopeOut, intercept;
};
#endif
