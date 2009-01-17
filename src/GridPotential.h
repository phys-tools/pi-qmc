// $Id: GridPotential.h,v 1.6 2008/11/23 21:29:50 jshumwa Exp $
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
#ifndef __GridPotential_h 
#define __GridPotential_h
class MultiLevelSampler;
class Paths;
class SimulationInfo;
#include "Action.h"
#include <string>
#include <vector>
#include <blitz/array.h>

/** Class for getting action and potential from a grid.
  * @todo Set up a potential table.
  * @bug Hard coded for NDIM=3.
  * @version $Revision: 1.6 $
  * @author John Shumway. */
class GridPotential : public Action {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,3> Array3;
  /// Constructor by providing an HDF5 file name and SimulationInfo.
  GridPotential(const SimulationInfo& simInfo, const std::string& filename);
  /// Virtual destructor.
  virtual ~GridPotential() {}
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
          double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
  /// The timestep.
  const double tau;
  /// The grid dimensions.
  int n1,n2,n3;
  /// Inverse grid spacing.
  double b;
  /// The potential grids.
  Array3 vegrid, vhgrid;
  /// Evaluate the potential at a point.
  double v(Vec, const int ipart) const;
  /// The grid index that associates the particles and grids.
  std::vector <Array3*> vindex;
};
#endif
