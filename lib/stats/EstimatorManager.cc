// $Id: EstimatorManager.cc,v 1.5 2008/11/19 21:26:23 jshumwa Exp $
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "EstimatorManager.h"
#include "AsciiReportBuilder.h"
#include "H5ReportBuilder.h"
#include "MPIManager.h"
#include "StdoutReportBuilder.h"
#include <algorithm>
#include <functional>
#include <iostream>

EstimatorManager::EstimatorManager(const std::string& filename,
    MPIManager *mpi, const SimInfoWriter *simInfoWriter)
  : filename(filename), mpi(mpi), simInfoWriter(simInfoWriter) {
  int rank=0;
#ifdef ENABLE_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  if (rank==0) {
    builders.push_back(new H5ReportBuilder(filename,simInfoWriter));
    builders.push_back(new StdoutReportBuilder());
    builders.push_back(new AsciiReportBuilder("pimc.dat"));
  }
}

EstimatorManager::~EstimatorManager() {
  delete simInfoWriter;
  for (EstimatorIter i=estimator.begin(); i!=estimator.end(); ++i) delete *i;
  for (BuilderIter i=builders.begin(); i!=builders.end(); ++i) delete *i;
}

void EstimatorManager::startWritingGroup(const int nstep, 
                                         const std::string& name) {
  this->nstep=nstep; istep=0;
  int rank=0;
#ifdef ENABLE_MPI
  rank = MPI::COMM_WORLD.Get_rank();
#endif
  if (rank==0) {
    for (BuilderIter builder=builders.begin();
         builder!=builders.end(); ++builder) {
      (*builder)->startWritingGroup(*this);
    }
  }
}

void EstimatorManager::writeStep() {
  // Calculate averages over clones.
  for (EstimatorIter est=estimator.begin(); est!=estimator.end(); ++est) {
    (*est)->averageOverClones(mpi);
  };
  // Write step data to using EstimatorReportBuilder objects.
  if (!mpi || mpi->isMain()) {
    for (BuilderIter builder=builders.begin();
         builder!=builders.end(); ++builder) {
      (*builder)->writeStep(*this);
    }
  }
  // Avoid a race condition.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Barrier();
    if (mpi->isCloneMain())  {
      mpi->getCloneComm().Barrier();
    }
  }
#endif
}

std::vector<Estimator*>& 
EstimatorManager::getEstimatorSet(const std::string& name) {
  std::vector<Estimator*>& set(estimatorSet[name]);
  if (name=="all" && set.size()!=estimator.size()) {
    set.resize(estimator.size());
    std::list<Estimator*>::iterator iter=estimator.begin();
    for (unsigned int i=0; i<set.size(); ++i) set[i]=*iter++;
  }
  return set;
}
