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
#ifndef __ArrayBlockedEstimator_h_
#define __ArrayBlockedEstimator_h_
#include "Estimator.h"
#include "EstimatorReportBuilder.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
class Paths;
/// Base class for scalar estimators.
/// @version $Revision$
/// @author John Shumway
class ArrayBlockedEstimator : public Estimator {
public:
  /// Construct by giving the name.
  ArrayBlockedEstimator(const std::string &name, 
                        const std::string &typeString, bool hasError);
  /// Virtual destructor.
  virtual ~ArrayBlockedEstimator();
  /// Clear the value.
  virtual void reset()=0;
  /// Callback EstimatorReportBuilder method for blocked arrays.
  virtual void startReport(EstimatorReportBuilder& builder) {
    builder.startArrayBlockedReport(*this);}
  /// Callback EstimatorReportBuilder method for blocked arrays.
  virtual void reportStep(EstimatorReportBuilder& builder) {
    builder.reportArrayBlockedStep(*this);}
  /// Get the number of dimensions in the array.
  virtual int getNDim() const=0;
  /// Get the extent of the array in dimensions idim.
  virtual int getExtent(const int idim) const=0;
  /// Get a pointer to the data.
  virtual const void* getData() const=0;
  virtual const void* getError() const=0;
  bool hasError() const {return hasErrorFlag;}
  virtual void normalize() const=0;
  virtual void unnormalize() const=0;
  virtual bool hasScale() const=0;
  virtual bool hasOrigin() const=0;
  virtual const void* getScale() const=0;
  virtual const void* getOrigin() const=0;
  virtual bool hasMin() const=0;
  virtual bool hasMax() const=0;
  virtual const void* getMin() const=0;
  virtual const void* getMax() const=0;
  protected:
  bool hasErrorFlag;
};
#endif
