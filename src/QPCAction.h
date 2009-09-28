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
#ifndef __QPCAction_h_
#define __QPCAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
template <int TDIM> class Beads;
#include <blitz/array.h>
#include "Action.h"

/** Class for calculating the action for a model of a quantum point contact.
  * From Hirose, Meir, and Wingreen PRL 90, 026804 (2002),
  * @f[ V_{\text{QPC}}(x,y)
  * =V(x)/2 + m^*(\omega_y+V(x)/\hbar]y^2/2 - V_{\text{wire}}(y),  @f]
  * where @f$ V(x)=V_0/\cosh^2(x/d) @f$
  * and @f$ V_{\text{wire}}(y)=(1/2)m^*\omega_y^2 y^2 @f$.
  *
  * Based on Hirose et al, we set the default values to
  * @f$ m^*=0.0667 @f$, @f$ \hbar\omega_y = 2.0 @f$ meV,
  * @f$ V_0 = 3.0 @f$ meV, and @f$ d=41.3 @f$ nm.
  *
  * @version $Revision$
  * @author John Shumway. */
class QPCAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor by providing the timestep tau.
  QPCAction(const double tau, 
            const double length=780.46, const double v0=0.00011025,
            const double omega=0.0000735, const double mass=0.0667);
  /// Virtual destructor.
  virtual ~QPCAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
 virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
private:
  /// The timestep.
  const double tau;
  /// The length of the QPC barrier.
  const double length;
  /// The height of the QPC barrier.
  const double v0;
  /// The lateral frequency of the quantum wire.
  const double omega;
  /// The mass.
  const double mass;
  /// Evaluate the potential at a point.
  double v(double x, double y) const;
};
#endif
