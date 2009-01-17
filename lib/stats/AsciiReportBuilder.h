// $Id: AsciiReportBuilder.h,v 1.2 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __AsciiReportBuilder_h_
#define __AsciiReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include <fstream>
#include <string>
#include <vector>

/** Class for reporting estimators to an ascii column file.
@version $Revision: 1.2 $
@author John Shumway */
class AsciiReportBuilder : public EstimatorReportBuilder {
public:
  /// Constructor.
  AsciiReportBuilder(const std::string & filename);
  /// Virtual destructor.
  virtual ~AsciiReportBuilder();
  /// Start a group for writing.
  virtual void startWritingGroup(EstimatorManager&);
  /// Write a simulation step to the current writing group.
  virtual void writeStep(EstimatorManager&);
  /// Method to write a step for a ScalarEstimator.
  virtual void reportScalarStep(const ScalarEstimator& est);
private:
  std::string filename;
  std::ofstream file;
  int nstep;
  int istep;
};
#endif
