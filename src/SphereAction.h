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
#ifndef __SphereAction_h_
#define __SphereAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
class Species;
template <int TDIM> class Beads;
#include "Action.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for calculating primitive action for a hard bounding sphere.
  * @todo May want to improve action for paths that go near the boundary.
  * @bug Only puts a sphere around the first particle of the species.
  * @todo Do a set OR on a set of spheres.
  * @version $Revision$
  * @author John Shumway. */
class SphereAction : public Action {
public:
  /// Typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<int,1> BArray;
  /// Constructor by providing the timestep tau.
  SphereAction(const double tau, const double radius, const Species&);
  /// Virtual destructor.
  virtual ~SphereAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
 virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int) const {return 0;}
  /// Calculate action and derivatives at a bead (no contribution).
  virtual void getBeadAction(const Paths&, int ipart, int islice,
      double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
    u=utau=ulambda=0; fm=0.; fp=0.;}
private:
  /// The timestep.
  const double tau;
  /// The radius.
  const double radius;
  /// The number of centers.
  const int ncenter;
  /// The center particles' indices.
  IArray ipart;
  /// Flags to see if center particles are moving.
  BArray isMoving;
  /// Indicies of moving centers.
  IArray icm;
};
#endif
