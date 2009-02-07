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
#ifndef __StdoutReportBuilder_h_
#define __StdoutReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include <string>
#include <vector>
#include <blitz/array.h>

/** Class for reporting estimators to standard output.
@version $Revision$
@author John Shumway */
class StdoutReportBuilder : public EstimatorReportBuilder {
public:
  /// Constructor.
  StdoutReportBuilder() {}
  /// Virtual destructor.
  virtual ~StdoutReportBuilder() {}
  /// Start a group for writing.
  virtual void startWritingGroup(EstimatorManager&);
  /// Write a simulation step to the current writing group.
  virtual void writeStep(EstimatorManager&);
  /// Report a step of scalar data.
  virtual void reportScalarStep(const ScalarEstimator& est);
  /// Method to start reporting a AccRejEstimator.
  virtual void startAccRejReport(const AccRejEstimator& est) {}
  /// Method to write a step for a AccRejEstimator.
  virtual void reportAccRejStep(const AccRejEstimator& est);
  /// Report a step of array data.
  virtual void reportArrayBlockedStep(const ArrayBlockedEstimator& est);
private:
  typedef blitz::Array<double,1> Array;
  int nstep;
  int istep;
  int iscalar;
  /// The sum of the estimator.
  Array sum, sum2, norm;
};
#endif
