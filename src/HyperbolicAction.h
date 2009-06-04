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
#ifndef __HyperbolicAction_h_
#define __HyperbolicAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
class PeriodicGaussian;
#include "Action.h"
#include <blitz/array.h>
#include <vector>

/** Class for calculating the free particle hyperbolic action.
 *
 * The effective mass increases with energy,
 * @f[ \epsilon(1+\alpha\epsilon) = \frac{k^2}{2m}, @f]
 * so the energy dispersion relation is hyperbolic, like in special relativity, 
 * @f[ \epsilon(k) = \frac{1}{2\alpha}\sqrt{1+\frac{2\alpha k^2}{m^*}}
 *     -\frac{1}{2\alpha}. @f] 
 * The propagator is
 * @f[ G(r,r';\tau) = \frac{e^{\frac{\tau}{2\alpha}}\tau}{32\pi^2\alpha z^2}
 *     \left(\frac{2m^*}{\alpha}\right)^{\frac{3}{2}} K_2(z),@f]
 * where @f$ K_2 @f$ is the modified Bessel function of the second kind
 * and 
 * @f[ z=\frac{\tau}{2\alpha}\sqrt{1+\frac{2m^*|r-r'|^2}{\tau^2}}. @f]
 *
 * @version $Revision$
 * @bug Periodic images are ignored.
 * @author John Shumway */
class HyperbolicAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<PeriodicGaussian*,3> PGArray;
  /// Construct by providing simulation info.
  HyperbolicAction(const SimulationInfo&, const int maxlevel);
  /// Virtual destructor.
  virtual ~HyperbolicAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
 virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead (defaults to no
  /// contribution).
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
};
#endif
