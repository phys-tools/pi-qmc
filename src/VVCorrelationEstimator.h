// $Id: VVCorrelationEstimator.h,v 1.10 2006/10/18 17:08:19 jshumwa Exp $
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
#ifndef __VVCorrelationEstimator_h_
#define __VVCorrelationEstimator_h_
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Action.h"
#include "Paths.h"
#include "SuperCell.h"

class Paths;
class Action;
class SimulationInfo;
class SuperCell;
/** @f$V \dot V@f$ Correlation Estimator.
 *  @version $Revision: 1.10 $
 *  @author Matthew Harowitz  */
class VVCorrelationEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VecArray;
  /// Constructor.
  VVCorrelationEstimator(const SimulationInfo& simInfo);
  /// Virtual destructor.
  virtual ~VVCorrelationEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start,
                          const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of momentum estimator.
  virtual double calcValue() {return value=vtot/vnorm;}
  /// Clear value of momentum estimator.
  virtual void reset() {vtot=vnorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// Array of velocities.
  VecArray v;
  /// @f$v(t) \dot v(t+\frac{\beta}{2})@f$
  double vdot;
  /// The accumulated velocity.
  double vtot;
  /// The normalization of the accumualted velocity.
  int vnorm;
  /// The supercell.
  SuperCell& cell;
};

#endif
