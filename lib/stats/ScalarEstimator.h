// $Id$
/*  Copyright (C) 2004-2007 John B. Shumway, Jr.

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
#ifndef __ScalarEstimator_h_
#define __ScalarEstimator_h_
#include "Estimator.h"
#include "EstimatorReportBuilder.h"
#include <blitz/tinyvec.h>
class Paths;
/// Base class for scalar estimators.
/// @version $Revision$
/// @author John Shumway
class ScalarEstimator : public Estimator {
public:
  /// Construct by giving the name.
  ScalarEstimator(const std::string& name);
  /// Construct by giving the name.
  ScalarEstimator(const std::string &name, const std::string &typeString,
                  const std::string &unitName, double shift, double scale);
  /// Virtual destructor.
  virtual ~ScalarEstimator() {}
  /// Get the value.
  virtual double calcValue()=0;
  /// Get the value.
  double getValue() const {return scale*(value+shift);}
  /// Get the value.
  void setValue(const double v) {value=v;}
  /// Clear the value.
  virtual void reset()=0;
  /// Average over clones.
  virtual void averageOverClones(const MPIManager* mpi);
  /// Callback EstimatorReportBuilder method for scalars.
  virtual void startReport(EstimatorReportBuilder& builder) {
    builder.startScalarReport(*this);}
  /// Callback EstimatorReportBuilder method for scalars.
  virtual void reportStep(EstimatorReportBuilder& builder) {
    builder.reportScalarStep(*this);}
protected:
  double value;
  /// Scale value.
  const double scale;
  /// Shift value.
  const double shift;
};
#endif
