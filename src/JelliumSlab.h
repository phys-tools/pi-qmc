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
#ifndef __JelliumSlab_h_
#define __JelliumSlab_h_
class SectionSamplerInterface;
class DisplaceMoveSampler;
template <int TDIM> class Beads;
class SimulationInfo;
#include "Action.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for calculating the primitive action for a jellium slab.
  * Can be used to run PIMC simulations like the DMC simulations
  * of <a href="http://link.aps.org/abstract/PRB/v54/p17199">
  * Acioli and Ceperley, PRB <b>54</b>, 17199-17207 (1996)</a>.
  * The slab is normal to the x-direction and has a thickness s.
  * For periodic boundary conditions, we need to include a uniform
  * background charge. The density of Q positive charges in the slab is
  * rhoSlab=Q/(s*Ly*Lz). The uniform background 
  * charge density is rhoBack=-Q/(Lx*Ly*Lz). The net density 
  * is thus rhoSlab+rhoBack inside the slab and rhoBack outside the slab.
  * The boolean flag isSlabInMiddle determines whether the slab or
  * barrier is in the middle of the cell. Call the respective densities
  * rhoMid and rhoEdge, and let w be the width of the middle region.
  * The solution to the Poisson equation is
  * @f[
  * \phi(x) = \begin{cases}
  * -2\pi \rho_{\text{mid}}x^2 + \phi_{\text{mid}}; & |x|<w/2 \\
  * -2\pi \rho_{\text{edge}}(L_x/2-|x|)^2 + \phi_{\text{edge}}; & |x|>w/2.
  * \end{cases}
  * @f]
  * The constants are chosen to make the potential continuous and 
  * sum to zero across the cell,
  * @f[
  * \begin{split}
  * \phi_{\text{mid}} &= \left(\frac{L_x\phi_{\text{sum}}}{L_x-w}
  *     +\phi_{\text{diff}}\right)\Big/\left(1+\frac{w}{L_x-w}\right)\\
  * \phi_{\text{edge}} &= \left(\frac{L_x\phi_{\text{sum}}}{w}
  *     -\phi_{\text{diff}}\right)\Big/\left(1+\frac{L_x-w}{w}\right)\\
  * \phi_{\text{sum}} &= \frac{4\pi}{3L_x}\left[n_{\text{mid}}
  *      \left(\frac{w}{2}\right)^2+n_{\text{edge}}
         \left(\frac{L_x-w}{2}\right)^2 \right]\\
  * \phi_{\text{diff}} &= 2\pi\left[ n_{\text{mid}}
         \left(\frac{w}{2}\right)^2 - n_{\text{edge}}
         \left(\frac{L_x-w}{2}\right)^2\right].
  * \end{split}
  * @f]
  * @version $Revision$
  * @author John Shumway. */
class JelliumSlab : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Construct by providing the timestep tau.
  JelliumSlab(const SimulationInfo &simInfo, const double qtot,
              const double slabWidth, const bool isSlabInMiddle);
  /// Virtual destructor.
  virtual ~JelliumSlab() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const SectionSamplerInterface&,
                                     const int level);
 virtual double getActionDifference(const DisplaceMoveSampler&,
				    const int nMoving){ return 0;};
   /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const;
private:
  /// Calculate the electrostatic potential at a point x.
  double phi(double x) const {
    if (fabs(x)<0.5*w) {
      return -2*PI*rhoMid*x*x + phiMid;
    } else {
      return -2*PI*rhoEdge*(0.25*lx*lx-lx*fabs(x)+x*x) + phiEdge;
    }
  }
  /// The timestep.
  const double tau;
  const double q;
  const double lx;
  double w;
  double rhoMid, rhoEdge;
  double phiMid, phiEdge;
  static const double PI;
};
#endif
