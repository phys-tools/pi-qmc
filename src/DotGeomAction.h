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
#ifndef __DotGeomAction_h_
#define __DotGeomAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
template <int TDIM> class Beads;
#include <cstdlib>
#include <blitz/array.h>
#include "Action.h"

/** Class for calculating the action for a nanostructure (quantum dot)
  * with a specified geometry.
  * @version $Revision$
  * @author John Shumway. */
class DotGeomAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Structure class for describing nanostructure geometry.
  class Structure {
    virtual double getPotentialDiff(const VArray&, const VArray&, int n)=0;
  };
  /// Constructor. 
  DotGeomAction(const double tau);
  /// Virtual destructor.
  virtual ~DotGeomAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
 virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
};
#endif
