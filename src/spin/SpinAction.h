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
#ifndef __SpinAction_h_
#define __SpinAction_h_
class Paths;
class MultiLevelSampler;
class DisplaceMoveSampler;
template <int TDIM> class Beads;
#include "../src/Action.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/array.h>
#include <vector>

/** Action class for spin action for B-field.
 *  @f[   \langle \uparrow | \sigma_{z} |\uparrow \rangle 
 *	 = \int d^{4}s \chi^{*}_{\downarrow}(\vec s) V(s) 
 * 	 \chi(\vec s)_{\uparrow}(s) = 1  @f]
 *  where @f[ V(\vec s)=6(1+2\alpha)^3 (\dfrac{s_{0}^{2}-s_{1}^{2}-s_{2}^{2}+s_{3}^{2}}{2})
 *                 (\dfrac{e^{-\alpha s^2}}{s^2} )@f]
 *  The Action for B-field( @f$ \vec {B}=B_0 \hat{z} @f$ ) is
 *  @f[ \dfrac{B_0}{2}\tau (V(s)+V(s')) @f]
 *    
 *
 *                 
 * @version $Revision$
 * @author John Shumway. */
class SpinAction : public Action {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,4> SVec;
  typedef blitz::TinyVector<double,3> Vec3;
  /// Constructor.
  SpinAction(const double tau, const double mass, const double bx,
     const double by, const double bz, const double omega, const double gc);
  /// Virtual destructor.
  virtual ~SpinAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  virtual double getActionDifference(const DisplaceMoveSampler&,
                                     const int nMoving);
   /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
       double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
  /// Accept last move.
  virtual void acceptLastMove();
private:
  /// The timestep.
  const double tau;
  /// The mass.
  const double mass;
  /// The magnetic field.
  const double bx,by,bz,gc,invgc;
  const double omega;
  const double l2;
};
#endif
