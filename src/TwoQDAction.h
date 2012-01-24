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
#ifndef __TwoQDAction_h_
#define __TwoQDAction_h_
class SectionSamplerInterface;class DisplaceMoveSampler;
template <int TDIM> class Beads;
#include "Action.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for calculating the action for two quantum dots at interdot 
 * distance d.
  * @f[ V(x) = \frac{1}{2}m \omega min((x-d/2)^2 +  y^2 ,
  *   (x+d/2)^2+y^2 )@f]
  * @f[ S(x,x') = \frac{k\tau}{4}({x'}^2+x^2) @f]
  * @version $Revision$
  * @author Daejin Shin. */
class TwoQDAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Construct by providing the timestep tau.
  TwoQDAction(const double tau, const double omega, const double mass,
              const double d, const double alpha);
  /// Virtual destructor.
  virtual ~TwoQDAction() {}
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
  /// The frequency.
  const double omega;
  /// The mass.
  const double mass;
  /// The iterdot distance.
  const double d; 
  /// The x-y anisotropy
  const double alpha;
};
#endif
