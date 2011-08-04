// $Id$
/*  Copyright (C) 2007 John B. Shumway, Jr.

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
#ifndef __StillWebAction_h_
#define __StillWebAction_h_
class MultiLevelSampler;class DisplaceMoveSampler;
template <int TDIM> class Beads;
#include "Action.h"
class SimulationInfo;
#include <cstdlib>
#include <blitz/array.h>
#include <vector>

/** Class for calculating the primitive action for the Stillinger-Weber
  * potential.
  *
  * Taken from Stillinger and Weber, PRB 31, 5262-5271 (1985),
  * online at http://link.aps.org/abstract/PRB/v31/p5262.
  * The pairwise interaction is a function of separation,
  * @f[
  * f_2(r)= A(Br^{-p}-r^{-q})\exp[(r-a)^{-1}], r<a
  * @f]
  * and zero for larger separations.
  *
  * There is an angular interaction between bonds,
  * @f[
  * h(r_{ij},r_{ik},\theta_{jik}) = \lambda\exp[\gamma(r_{ij}-a)^{-1} +
  * \gamma(r_{ik}-a)^{-1}](\cos\theta_{jik}+\frac{1}{3})^2.
  * @f] 
  *
  * Si-Si lattice paramters is 10.246 Bohr radii (T=0K). 
  *
  * @bug Assumes only one particle is moving.
  * @todo Add neighbor table.
  * @version $Revision$
  * @author John Shumway. */
class StillWebAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,2> Array2;
  /// Construct by providing the simulation info and structure
  /// file (for neighbor tables).
  StillWebAction(const SimulationInfo&, const std::string &filename);
  /// Virtual destructor.
  virtual ~StillWebAction() {}
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

  /// Evaluate the pair potential.
  virtual double f2(const double r, const int i, const int j) const;
  /// Evaluate the angle potential.
  virtual double h(const double r1, const double r2, const double costheta,
           const int i, const int j, const int k) const;


private:
  /// Subclass for storage of species dependent parameters.
  class SWParam {
  public:
    double lambda;
    double sigmainv;
    double epsilon;
  };
  /// The timestep.
  const double tau;
  /// The number of particles.
  const int npart;
  /// SW parameters.
  double A;
  double B;
  double p;
  double q;
  double a;
  double gamma;
  /// Parameters for SiSi, GeGe, SiGe. 
  std::vector<SWParam> param;
  /// Neighbor tables.
  Array2 nn;
  /// Maximum number of neighbors in neighbor table.
  int nnmax;
};
#endif
