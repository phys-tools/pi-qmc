// $Id$
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#ifndef __WellImageAction_h_
#define __WellImageAction_h_
class SectionSamplerInterface;
class DisplaceMoveSampler;
template <int TDIM> class Beads;
class SimulationInfo;
class SuperCell;
#include "Action.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for calculating the primitive action for image charges around
  * a dielectric slab. The image charge potential is described in 
  * <a href="http://link.aps.org/abstract/PRB/v40/p12359">M. Kumagai, 
  * and T. Takagahara, "Excitonic and non-linear optical properties of
  * dielectric quantum-well structures," PRB 40, 12359-12381, (1989)</a>.
  *
  * The potential is given by Eqs. (2.8)-(2.15) of that paper.
  * Let @f$ z_1 @f$ and @f$ z_2 @f$ be the coordinates of two particles
  * normal to the center of the well, and let @f$ \rho @f$ be the component
  * of their separation vector tangental to the plane of the well.
  * Also let L, C, R, refer to the region where the two particles reside.
  * Define the following functions of the dielctrics,
  * @f[
  * p_1 = \frac{\epsilon_{in}-\epsilon_{out}}{\epsilon_{in}+\epsilon_{out}},
  * @f]
  * and 
  * @f[
  * p_2 = \frac{2\epsilon_{in}}{\epsilon_{in}+\epsilon_{out}}.
  * @f]
  * Each particle feels the effect of an image charge, or self energy,
  * @f[
  * V_{self}^C(z) = \frac{p_1 q^2}{2\epsilon_{in}|2z + L + 2\delta|}
  *                +\frac{p_1 q^2}{2\epsilon_{in}|2z - L - 2\delta|}
  *                +\sum_{n=2}^\infty \frac{p_1^n q^2}{2\epsilon_{in}}
  *                        \left[\frac{1}{|z - (-1)^n z + nL|}
  *                             +\frac{1}{|z - (-1)^n z - nL|}\right],
  * @f]
  * @f[
  * V_{self}^R(z) = V_{self}^L(-z)
  *               = -\frac{p_1 q^2}{2\epsilon_{out}|2z - L + 2\delta|}
  *                +p_2\sum_{n\ odd}\frac{p_1^n q^2}
  *                            {(\epsilon_{in}+\epsilon_{out})|2z+nL|}
  * @f]
  * The Coulomb interactions become
  * @f[
  * V^{CC}(z_1,z_2,\rho) =\frac{q_1q_2}{\epsilon_{in}\sqrt{\rho^2+(z_1-z_2)^2}}
  *   +\sum_{n=1}^\infty\frac{p_1^n q_1 q_2}{\epsilon_{in}}
  *          \left[\frac{1}{\sqrt{\rho^2+[z_1-(-1)^n z_2 -nL]^2}}
  *               +\frac{1}{\sqrt{\rho^2+[z_1-(-1)^n z_2 +nL]^2}}\right],
  * @f]
  * @f[
  * V^{LC}(z_1,z_2,\rho) = p_2\sum_{n=0}^\infty
  *     \frac{p_1^n q_1 q_2}{\epsilon_{in}\sqrt{\rho^2+[z_1-(-1)^nz_2-nL]^2}},
  * @f]
  * @f[
  * V^{RC}(z_1,z_2,\rho) = p_2\sum_{n=0}^\infty
  *     \frac{p_1^n q_1 q_2}{\epsilon_{in}\sqrt{\rho^2+[z_1-(-1)^nz_2+nL]^2}},
  * @f]
  * @f[
  * V^{LR}(z_1,z_2,\rho) = p_2^2\sum_{n=0}^\infty
  *     \frac{p_1^{2n} q_1 q_2}{\epsilon_{in}\sqrt{\rho^2+[z_1-z_2+2nL]^2}},
  * @f]
  * @f[
  * V^{RR}(z_1,z_2,\rho) =\frac{q_1q_2}{\epsilon_{out}\sqrt{\rho^2+(z_1-z_2)^2}}
  *          -\frac{p_1q_1q_2}{\epsilon_{out}\sqrt{\rho^2+(z_1+z_2-L)^2}}
  *  +p_2^2\sum_{n\ odd}
  *     \frac{p_1^{n} q_1 q_2}{\epsilon_{in}\sqrt{\rho^2+[z_1+z_2+nL]^2}}.
  * @f]
  * We will omit the direct terms (that diverge) since they are already 
  * handled in CoulombAction. This is an approximation, since CoulombAction
  * does not know about the different dielectric constants.
  * @version $Revision$
  * @author John Shumway. */
class WellImageAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  /// Construct by providing the timestep tau.
  WellImageAction(const SimulationInfo &simInfo, double epsIn,
    double epsOut, double width, double z0, double delta);
  /// Virtual destructor.
  virtual ~WellImageAction() {}
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
  /// Simulation timestep.
  const double tau;
  /// Number of particles.
  const int npart;
  /// Dielectric constants.
  const double epsIn, epsOut;
  /// Well width.
  const double width;
  /// Well center coordinate.
  const double z0;
  /// Distance delta for regularizing self energy divergence.
  const double delta;
  /// Dielectric factors.
  const double p1, p2;
  /// Particle charges.
  Array q;
  /// Self energy function.
  inline double vSelfC(double z) const;
  /// Self energy function.
  inline double vSelfR(double z) const;
  /// Coulomb interaction function.
  inline double vRC(double z1, double z2, double rho) const;
  /// Coulomb interaction function.
  inline double vRR(double z1, double z2, double rho) const;
  /// Coulomb interaction function.
  inline double vLR(double z1, double z2, double rho) const;
  /// Coulomb interaction function.
  inline double vCC(double z1, double z2, double rho) const;
  /// The SuperCell.
  const SuperCell &cell;
};
#endif
