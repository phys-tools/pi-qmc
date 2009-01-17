// $Id: SpringAction.h,v 1.9 2007/10/03 12:53:56 jshumwa Exp $
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
#ifndef __SpringAction_h_
#define __SpringAction_h_
class MultiLevelSampler;
class Paths;
class SimulationInfo;
class PeriodicGaussian;
#include "Action.h"
#include <blitz/array.h>
#include <vector>

/** Class for calculating the free particle kinetic ``spring'' action.
 * 
 * The kinetic action for a single particle is 
 * @f$ u(r,r';\tau) = \frac{Nd}{2}\log 2\pi\tau/m + m|r-r'|^2/(2\tau) @f$.
 *
 * @version $Revision: 1.9 $
 * @author John Shumway */
class SpringAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<bool,1> BArray;
  typedef blitz::Array<PeriodicGaussian*,3> PGArray;
  /// Construct by providing simulation info.
  SpringAction(const SimulationInfo&, const int maxlevel, const double pgDelta);
  /// Virtual destructor.
  virtual ~SpringAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead (defaults to no
  /// contribution).
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The inverse particle mass, @f$\lambda=1/2m@f$.
  Array lambda;
  /// The timestep.
  const double tau;
  PGArray pg;
  IArray specIndex;
  /// Flag for static particles.
  BArray isStatic;
};
#endif
