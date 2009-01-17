// $Id: H5ReportBuilder.h,v 1.5 2008/11/21 22:19:12 jshumwa Exp $
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
#ifndef __H5ReportBuilder_h_
#define __H5ReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include <string>
#include <vector>
#include <H5Cpp.h>
#include <hdf5.h>

#include "EstimatorManager.h"

/** Class for recording statistical data to an HDF5 file.

The H5ReportBuilder s responsible for recording data
recording it to disk. We have choosen an HDF5
file format (http://hdf.ncsa.uiuc.edu/HDF5) so that we can compress
the data and include meta data and other structure.
@version $Revision: 1.5 $
@author John Shumway */
class H5ReportBuilder : public EstimatorReportBuilder {
public:
  /// Constructor.
  H5ReportBuilder(const std::string &filename, 
                  const EstimatorManager::SimInfoWriter*);
  /// Virtual destructor.
  virtual ~H5ReportBuilder();
  /// Start a group for writing.
  virtual void startWritingGroup(EstimatorManager&);
  /// Write a simulation step to the current writing group.
  virtual void writeStep(EstimatorManager&);
  /// Method to start reporting a ScalarEstimator.
  virtual void startScalarReport(const ScalarEstimator& est);
  /// Method to write a step for a ScalarEstimator.
  virtual void reportScalarStep(const ScalarEstimator& est);
  /// Method to start reporting a ArrayBlockedEstimator.
  virtual void startArrayBlockedReport(const ArrayBlockedEstimator& est);
  /// Method to write a step for a ArrayBlockedEstimator.
  virtual void reportArrayBlockedStep(const ArrayBlockedEstimator& est);
  /// Method to start reporting a AccRejEstimator.
  virtual void startAccRejReport(const AccRejEstimator& est);
  /// Method to write a step for a AccRejEstimator.
  virtual void reportAccRejStep(const AccRejEstimator& est);
private:
  std::string filename;
  const EstimatorManager::SimInfoWriter *simInfoWriter;
  int nstep;
  int istep;
  typedef std::vector<hid_t> DataSetContainer;
  //typedef std::vector<H5::DataSet*> DataSetContainer;
  typedef DataSetContainer::iterator DataSetIter;
  /// The HDF5 file. 
  hid_t fileID;
  //H5::H5File file;
  /// The current HDF5 writing group.
  hid_t writingGroupID;
  //H5::Group* writingGroup;
  /// The group's step counter.
  hid_t stepAttrID;
  //H5::Attribute* stepAttr;
  /// The current HDF5 datasets.
  DataSetContainer dataset;
  /// Iterator pointing at the current dataset.
  DataSetIter dset;
};
#endif
