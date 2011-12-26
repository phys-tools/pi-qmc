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
#ifndef __ArrayEstimator_h_
#define __ArrayEstimator_h_
#include "Estimator.h"
#include "EstimatorReportBuilder.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
class Paths;
/// Base class for scalar estimators.
/// @version $Revision$
/// @author John Shumway
class ArrayEstimator : public Estimator {
public:
  /// Construct by giving the name.
  ArrayEstimator(const std::string& name) : Estimator(name) {}
  /// Virtual destructor.
  virtual ~ArrayEstimator() {}
  /// Get the value.
  virtual double calcValue()=0;
  /// Clear the value.
  virtual void reset()=0;
  /// Average over clones.
  virtual void averageOverClones(const MPIManager* mpi)=0;
  /// Callback EstimatorReportBuilder method for scalars.
  virtual void startReport(EstimatorReportBuilder& builder) {
    builder.startArrayReport(*this);}
  /// Callback EstimatorReportBuilder method for scalars.
  virtual void reportStep(EstimatorReportBuilder& builder) {
   builder.reportArrayStep(*this);}
  /// Get the number of dimensions in the array.
  virtual int getNDim() const=0;
  /// Get the extent of the array in dimensions idim.
  virtual int getExtent(const int idim) const=0;
protected:
};
#endif
