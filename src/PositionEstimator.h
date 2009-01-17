//$Id:
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

#ifndef __PositionEstimator_h_
#define __PositionEstimator_h_
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <iostream>
class Paths;
class SimulationInfo;
/** position estimator.
 *  @author Matthew Harowitz  */
class PositionEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<std::string,1> SArray;
  /// Constructor.
  PositionEstimator(const SimulationInfo& simInfo, const Species& species, const int index=2);
  /// Virtual destructor.
  virtual ~PositionEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of position estimator.
  virtual double calcValue() {return value=ptot/pnorm;}
  /// Clear value of position estimator.
  virtual void reset() {ptot=pnorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// Position .
  double position;
  /// The accumulated polarization.
  double ptot;
  /// The normalization of the accumualted polarization.
  double pnorm;
  /// The component of the position to measure (0->x, 1->y, 2->z)
  int index;
  /// The names of the particles.
  SArray names;
  /// The species to measure
  std::string spec;
};

#endif
