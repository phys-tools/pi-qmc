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
#ifndef __ProbDensityGrid_h_
#define __ProbDensityGrid_h_

#include <vector>
#include "Algorithm.h"
#include "LinkSummable.h"
#include <blitz/array.h>
class SimulationInfo;
class Paths;

/// Class for storing probability density on grids.
/// @version $Revision$
/// @author John Shumway 
class ProbDensityGrid : public Algorithm, public LinkSummable {
public:
  /// Constants and typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::Array<long int,NDIM> LIArray;
  /// Constructor.
  ProbDensityGrid(const IVec n, const double a,
                  const SimulationInfo&, const Paths*);
  /// Virtual destructor.
  virtual ~ProbDensityGrid();
  /// Run method.
  virtual void run() {;}
  /// Get the grid for species i.
  LIArray& getGrid(const int i) const {return *grid[i];}
  /// Get the current normalization.
  double getNorm() const {return norm;}
  /// Get the lattice constant.
  double getLatticeConstant() const {return a;}
  /// Get the grid dimensions.
  IVec getGridDim() const {return n;}
  /// Add the current path positions to bins.
  void bin();
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
protected:
  /// The grids for storing probability densities.
  std::vector<LIArray*> grid;
  /// The path storage.
  const Paths* paths;
  /// The normalization. 
  double norm;
  /// Dimensions of the grids.
  const IVec n;
  /// The lattice spacing.
  const double a;
  /// Reciprical lattice vectors.
  std::vector<Vec> b;
  /// Index from particle to species.
  std::vector<int> index;
};
#endif
