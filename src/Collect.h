// $Id: Collect.h,v 1.6 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __Collect_h_
#define __Collect_h_

class Estimator;
class EstimatorManager;
#include "Algorithm.h"
#include <vector>
#include <fstream>

/// Algorithm class for collecting estimator measurements.
/// @version $Revision: 1.6 $
/// @author John Shumway
class Collect : public Algorithm {
public:
  /// Constructor.
  Collect(std::string& estName, EstimatorManager&,  const int nstep);
  /// Virtual destructor.
  virtual ~Collect() {}
  /// Run the algorithm.
  virtual void run();
private:
  /// Reference to the EstimatorManager.
  EstimatorManager& estManager;
  /// Number of times this Collect object will be run.
  const int nstep;
  /// Internal flag for catching the first call to run.
  bool isFirstRun;
};
#endif
