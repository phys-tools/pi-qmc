// $Id$
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#ifndef __WindingEstimator_h_
#define __WindingEstimator_h_
#include "LinkSummable.h"
#include "stats/BlitzArrayBlkdEst.h"
#include <blitz/array.h>
#include "Species.h"

class Paths;
class SuperCell;
class SimulationInfo;
class MPIManager;
/** Estimator for winding around periodic boundary conditions.
 *  @version $Revision$
 *  @author John Shumway  */
class WindingEstimator : public BlitzArrayBlkdEst<NDIM>, public LinkSummable {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  WindingEstimator(const SimulationInfo& simInfo, int nmax,
    const std::string &name, bool isChargeCoupled, const Species *species, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~WindingEstimator();
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths);
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(int nslice);
private:
  ///
  int nmax;
  int npart;
  int ifirst;
  Vec winding;
  SuperCell &cell;
  IArray charge;
  bool isChargeCoupled;
  MPIManager *mpi;
};

#endif
