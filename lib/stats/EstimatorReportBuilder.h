// $Id: EstimatorReportBuilder.h,v 1.3 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __EstimatorReportBuilder_h_
#define __EstimatorReportBuilder_h_

class EstimatorManager;
class ScalarEstimator;
class AccRejEstimator;
class ArrayEstimator;
class ArrayBlockedEstimator;

/** Base class for reporting estimators.

We need to write methods for reporting each type of 
Estimator data.

@version $Revision: 1.3 $
@author John Shumway */
class EstimatorReportBuilder {
public:
  /// Constructor.
  EstimatorReportBuilder() {}
  /// Virtual destructor.
  virtual ~EstimatorReportBuilder() {}
  /// Start a group for writing.
  virtual void startWritingGroup(EstimatorManager&)=0;
  /// Write a simulation step to the current writing group.
  virtual void writeStep(EstimatorManager&)=0;
  /// Method to start reporting a ScalarEstimator.
  virtual void startScalarReport(const ScalarEstimator& est) {}
  /// Method to write a step for a ScalarEstimator.
  virtual void reportScalarStep(const ScalarEstimator& est) {}
  /// Method to start reporting a AccRejEstimator.
  virtual void startAccRejReport(const AccRejEstimator& est) {}
  /// Method to write a step for a AccRejEstimator.
  virtual void reportAccRejStep(const AccRejEstimator& est) {}
  /// Method to start reporting a ArrayEstimator.
  virtual void startArrayReport(const ArrayEstimator& est) {}
  /// Method to write a step for a ArrayEstimator.
  virtual void reportArrayStep(const ArrayEstimator& est) {}
  /// Method to start reporting a ArrayBlockedEstimator.
  virtual void startArrayBlockedReport(const ArrayBlockedEstimator& est) {}
  /// Method to write a step for a ArrayBlockedEstimator.
  virtual void reportArrayBlockedStep(const ArrayBlockedEstimator& est) {}
};
#endif
