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
#ifndef __SHODotAction_h_
#define __SHODotAction_h_
class SectionSamplerInterface;class DisplaceMoveSampler;
class Species;
template <int TDIM> class Beads;
#include <cstdlib>
#include <blitz/array.h>
#include "Action.h"

/** Class for calculating the action for a simple harmonic oscillator.
  * @version $Revision$
  * @author John Shumway. */
class SHODotAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor by providing the timestep tau, thickness t, potenial
  /// v0 and spring constant k.
  SHODotAction(double tau, double t, double v0, double omega, 
    double z, const Species&);
  /// Virtual destructor.
  virtual ~SHODotAction() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const SectionSamplerInterface&,
                                     int level);
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int iLastSlice);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
  /// The thickness.
  const double t;
  /// The well depth.
  const double v0;
  /// The spring constant.
  const double k;
  const double z;
  /// The first particle and particle count for this species.
  const int ifirst, npart;
};
#endif
