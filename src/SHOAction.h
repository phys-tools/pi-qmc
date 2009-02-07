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
#ifndef __SHOAction_h_
#define __SHOAction_h_
class MultiLevelSampler;
template <int TDIM> class Beads;
#include <blitz/array.h>
#include "Action.h"

/** Class for calculating the action for a simple harmonic oscillator.
  * The thermal density matrix for a simple harmonic oscillator is,
  * @f[ \rho(x,x';\tau) = \sqrt{\frac{m\omega}{2\pi\sinh\omega\tau}}
  *        \exp\frac{-m\omega[(x^2+{x'}^2)\cosh\omega\tau-2xx']
  *                   }{2\sinh\omega\tau}. @f]
  * For a derivation, see Feynman, <em>Statistical Mecahnics</em> (1972).
  * This class only represents the potential contribution to the
  * action, defined as @f$U=-\log(\rho/\rho^0)@f$.
  * This action is
  * @f[ U(x,x';\tau) = \frac{1}{2}\log\frac{2\pi\sinh\omega\tau}{m\omega}
  *  +\frac{m\omega[(x^2+{x'}^2)\cosh\omega\tau-2xx']}{2\sinh\omega\tau}
  *  -\frac{1}{2}\log 2\pi\tau/m -\frac{m(x-x')^2}{2\tau}. @f]
  * We also need the time derivative,
  * @f[ \dot{U}(x,x';\tau) =
  *   \frac{\omega}{2}\coth\omega\tau
  *  +\frac{m\omega^2}{2}\left[(x^2+{x'}^2)(1-\coth^2\omega\tau)
  *  +\frac{2xx'\cosh\omega\tau}{\sinh^2\omega\tau}\right]
  *  +\frac{m(x-x')^2}{2\tau^2}-\frac{1}{2\tau}. @f]
  *
  * For the virial estimator, we neee to calculate forces.
  * @f[ \nabla_{x}U(x,x';\tau)=
  *   \frac{m\omega[x\cosh\omega\tau-x']}{\sinh\omega\tau}
  *   -\frac{m(x-x')}{\tau} @f] 
  * (These formulas generalize to several dimensions by using vector
  * notation.)
  * @todo Treat masses correctly.
  * @todo Calculate forces for virial estimator.
  * @bug Doesn't use mass correctly.
  * @version $Revision$
  * @author John Shumway. */
class SHOAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor by providing the timestep tau.
  SHOAction(const double tau, const double omega, const double mass,
            const int ndim=NDIM);
  /// Virtual destructor.
  virtual ~SHOAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
  /// The frequency.
  const double omega;
  /// The mass.
  const double mass;
  /// The number of dimensions (can be less than NDIM, i.e., to make a wire).
  const int ndim;
};
#endif
