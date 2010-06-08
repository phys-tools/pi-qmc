// $Id: EMARateAction.h 185 2009-10-13 06:00:07Z john.shumwayjr $
/*  Copyright (C) 2010 John B. Shumway, Jr.

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
#ifndef __EMARateAction_h_
#define __EMARateAction_h_
class MultiLevelSampler;
class DisplaceMoveSampler;
class Paths;
class SimulationInfo;
class Species;
#include "Action.h"
#include <blitz/array.h>
#include <vector>

/**
 * @version $Revision: 185 $
 * @author John Shumway */
class EMARateAction : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<bool,1> BArray;
  /// Construct by providing simulation info.
  EMARateAction(const SimulationInfo&, const Species&, const Species&);
  /// Virtual destructor.
  virtual ~EMARateAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the difference in action.
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int nslice);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate action and derivatives at a bead (defaults to no
  /// contribution).
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
  /// The electron species.
  const Species& species1;
  /// The hole species.
  const Species& species2;
};
#endif
