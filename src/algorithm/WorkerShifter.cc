#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "WorkerShifter.h"
#include "base/Paths.h"
#include "stats/MPIManager.h"
#include "util/RandomNumGenerator.h"

WorkerShifter::WorkerShifter(const int maxShift, Paths& paths, MPIManager* mpi)
  : CompositeAlgorithm(1), maxShift(maxShift), paths(paths), mpi(mpi) {
}

WorkerShifter::~WorkerShifter() {}

void WorkerShifter::run() {
  int ishift=(int)((maxShift-1)*RandomNumGenerator::getRand()*(1-1e-8))+1;
  if (maxShift==0) ishift=0;
#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(&ishift,1,MPI::INT,0);
#endif
  paths.shift(ishift);
  CompositeAlgorithm::run();
}
