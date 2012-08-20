#include "config.h"
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
#include <fstream>

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
  //char buffer[256];
  //int namelength;
  //MPI::Get_processor_name (buffer, namelength);
  //std::string name(buffer, namelength); 
  //std::time_t rawtime;
  //std::time(&rawtime);
  //std::cout << "Barrier reached by " << name << " at " 
  //          << ctime(&rawtime) << std::endl;
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


void EstimatorManager::recordInputDocument(const std::string &filename) {
  std::string buffer;
  std::string docstring;
  std::ifstream in(filename.c_str());
  while(!std::getline(in, buffer).eof())
    docstring+=buffer+="\n";

  for (BuilderIter builder=builders.begin();
       builder!=builders.end(); ++builder) {
    (*builder)->recordInputDocument(docstring);
  }
}
