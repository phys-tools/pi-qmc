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
#ifndef __EstimatorManager_h_
#define __EstimatorManager_h_
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <hdf5.h>
#include "Estimator.h"
class MPIManager;
class EstimatorReportBuilder;

/** Class for managing estimators. 
Stores all estimators and sets of estimators.

We define an Estimator to be any of the mathematical expressions
that are sampled during the simulation. These can be physical
quantities, such as ThermalEnergyEstimator or ColoumbEnergyEstimator.
There are also Estimator classes for algorithmic quantities, such
accept/reject ratios.

The EstimatorManager is responsible for coordinating all the data
collection and recording it to disk. We have choosen an HDF5
file format (http://hdf.ncsa.uiuc.edu/HDF5) so that we can compress
the data and include meta data and other structure.

@todo Need to make public interface for access by EstimatorReportBuilder.
@version $Revision$
@author John Shumway */
class EstimatorManager{
public:
  /// Optional helper class for writting simulation info.
  class SimInfoWriter{
    public: 
    virtual void writeH5(hid_t) const {};
    virtual void writeStdout(std::ostream&) const {};
    virtual void writeAscii(std::ostream&) const {};
  };
  /// Constructor.
  EstimatorManager(const std::string& filename, MPIManager *mpi,
                   const SimInfoWriter*);
  /// Set estimator array to the named estimator set.
  std::vector<Estimator*>& getEstimatorSet(const std::string &name);
  /// Virtual destructor.
  virtual ~EstimatorManager();
  /// Create a new estimator set.
  void createEstimatorSet(const std::string &name, 
                          const std::vector<const std::string>&);
  /// Add an estimator.
  void add (Estimator* e) {estimator.push_back(e);}
  /// Start a group for writing.
  void startWritingGroup(const int nstep, const std::string& name);
  /// Write a simulation step to the current writing group.
  void writeStep();
  /// Record the input file in the simulation output.
  void recordInputDocument(const std::string &filename);
friend class H5ReportBuilder;
friend class StdoutReportBuilder;
friend class AsciiReportBuilder;
private:
  typedef std::list<Estimator*> EstimatorList;
  typedef EstimatorList::iterator EstimatorIter;
  typedef std::list<EstimatorReportBuilder*> BuilderList;
  typedef BuilderList::iterator BuilderIter;
  /// HDF5 output file name.
  const std::string filename;
  /// Estimators.
  EstimatorList estimator;
  /// Estimator sets.
  std::map<std::string,std::vector<Estimator*> > estimatorSet;
  /// The number of steps in the current writing group.
  int nstep;
  /// The current step number in the current writing group.
  int istep;
  /// The MPIManager.
  MPIManager *mpi;
  /// A pointer to the simulation info writer.
  const SimInfoWriter *simInfoWriter;
  /// The report builders.
  BuilderList builders;
};
#endif
