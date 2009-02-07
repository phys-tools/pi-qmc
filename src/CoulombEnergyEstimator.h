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
#ifndef __CoulombEnergyEstimator_h_
#define __CoulombEnergyEstimator_h_
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <blitz/array.h>
class Paths;
class Action;
class SimulationInfo;
class MPIManager;
/** Coulomb energy estimator. 
 *  @version $Revision$
 *  @author John Shumway  */
class CoulombEnergyEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor.
  CoulombEnergyEstimator(const SimulationInfo& simInfo, const Action*, 
       const double epsilon, MPIManager *mpi,
       const std::string& unitName, double scale, double shift);
  /// Virtual destructor.
  virtual ~CoulombEnergyEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of coulomb energy estimator.
  virtual double calcValue() {return value=etot/enorm;}
  /// Clear value of coulomb energy estimator.
  virtual void reset() {etot=enorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// CoulombEnergy energy.
  double energy;
  /// The accumulated energy.
  double etot;
  /// The normalization of the accumualted energy.
  double enorm;
  /// Action.
  const Action* action;
  ///dielectric constant
  const double epsilon;
  /// The charges of the particles.
  Array q;
  /// The MPI manager.
  const MPIManager* mpi;
};

#endif
