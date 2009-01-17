// $Id: Estimator.h,v 1.3 2007/10/26 07:04:26 jshumwa Exp $
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
#ifndef _Estimator_h__
#define _Estimator_h__

#include <string>
#include <iostream>
class MPIManager;
class Paths;
class EstimatorReportBuilder;
/** Base class for all estimators.

Estimators should be subclassed to handle different
data types.

@version $Revision: 1.3 $
@author John Shumway */
class Estimator {
public:
  /// Construct by giving the name of the estimator.
  Estimator(const std::string& name) : name(name) {}
  /// Construct by giving the name of the estimator.
  Estimator(const std::string& name, const std::string& unitName) 
    : name(name), unitName(unitName) {}
  /// Virtual destructor.
  virtual ~Estimator() {}
  /// Get the name of the estimator.
  const std::string& getName() const {return name;}
  /// Get the name of the unit for the estimator.
  const std::string& getUnitName() const {return unitName;}
  /// Evaluate estimator for a path configuration.
  virtual void evaluate(const Paths&)=0;
  /// Average over clones.
  virtual void averageOverClones(const MPIManager*) {};
  /// Callback EstimatorReportBuilder method.
  virtual void startReport(EstimatorReportBuilder& builder) {}
  /// Callback EstimatorReportBuilder method.
  virtual void reportStep(EstimatorReportBuilder& builder) {}
private:
  /// The name of the estimator.
  const std::string name;
  /// The unit name;
  const std::string unitName;
};
#endif
