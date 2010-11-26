// $Id: FreeEnergyEstimator.h 181 2009-09-28 22:48:51Z john.shumwayjr $
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
#ifndef __FreeEnergyEstimator_h_
#define __FreeEnergyEstimator_h_
#include "stats/BlitzArrayBlkdEst.h"
#include "Paths.h"
#include <fftw3.h>
#include <blitz/array.h>
#include "Species.h"
class Paths;
class SimulationInfo;
class MPIManager;
/** Permutation estimator for a homogeneous system.
 *  @version $Revision: 181 $
 *  @bug Hard coded for two particles and no MPI.
 *  @author John Shumway  */

class FreeEnergyEstimator : public BlitzArrayBlkdEst<1> {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<bool,1> BArray;
  typedef blitz::Array<double,1> Array;
  /// Constructor.
  FreeEnergyEstimator(const SimulationInfo& simInfo, int nmodel,
      MPIManager *mpi);
  /// Virtual destructor.
  virtual ~FreeEnergyEstimator();
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths);

  virtual void averageOverClones(const MPIManager* mpi);

private:
  ///
  MPIManager *mpi;
};

#endif
