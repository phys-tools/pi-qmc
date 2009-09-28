// $Id:
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
#ifndef __EFieldAction_h
#define __EFieldAction_h
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
class SuperCell;
#include "Action.h"
#include <blitz/array.h>

/** Class for getting action from stepwise electric field.
  * @bug Hard coded for NDIM=3.
  * @version $Revision:
  * @author Matthew Harowitz. */
class EFieldAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<double, 1> Array;
  typedef blitz::Array<int, 1> IArray;
  typedef blitz::TinyVector<double, NDIM> Vec;
  /// Constructor by providing SimulationInfo and a scale factor.
  EFieldAction(const SimulationInfo& simInfo, const double scale=0., const int index=2);
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
  /// The scale factor
  const double scale;
  /// The component of the electric field (0->x, 1->y, 2->z)
  const int component;
  /// The supercell
  SuperCell& cell;
  /// The electric potential
  double v(const double z) const;
  /// The size of the physical region.
  double a;
};
#endif
