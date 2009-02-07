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
#ifndef __WriteProbDensity_h_
#define __WriteProbDensity_h_

class ProbDensityGrid;
class SimulationInfo;
#include "Algorithm.h"
#include <string>
#include <vector>
#include <blitz/array.h>

/// Class for writing the probability density to a hdf5 file.
/// @version $Revision$
/// @author John Shumway
class WriteProbDensity : public Algorithm {
public:
  /// Constants and typedefs.
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<long int,NDIM> LIArrayN;
  typedef blitz::Array<float,NDIM> FArrayN;
  /// Constructor.
  WriteProbDensity(const SimulationInfo&, const ProbDensityGrid*, 
                   const std::string& filename);
  /// Virtual destructor.
  virtual ~WriteProbDensity();
  /// Algorithm run method.
  virtual void run();
private:
  /// Name of output file.
  const std::string filename;
  /// The grids.
  const ProbDensityGrid* grids;
  /// The names of the grids.
  std::vector<std::string> name;
};
#endif
