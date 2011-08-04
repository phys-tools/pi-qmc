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
#ifndef __ConditionalDensityGrid_h_
#define __ConditionalDensityGrid_h_

#include <vector>
#include "ProbDensityGrid.h"
#include <cstdlib>
#include <blitz/array.h>
class SimulationInfo;
class Paths;
class Species;

/// Class for storing conditional probability density on grids.
/// This is one way to get a correlation function in a three-dimensional
/// structure: only collect density when a particular particle is within
/// some radius of a specified point.
/// @version $Revision$
/// @author John Shumway 
class ConditionalDensityGrid : public ProbDensityGrid {
public:
  /// Constructor.
  ConditionalDensityGrid(const IVec n,
		  const double a, const SimulationInfo&, const Paths*,
		  const Species&, Vec center, const double radius);
  /// Virtual destructor.
  virtual ~ConditionalDensityGrid() {;}
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
private:
  /// Index of conditional particle.
  const int ifirst, ilast;
  /// Center of conditional region.
  const Vec center;
  /// Radius of the conditoinal region.
  const double radius;
};
#endif
