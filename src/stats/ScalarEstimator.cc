#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ScalarEstimator.h"
#include "MPIManager.h"

ScalarEstimator::ScalarEstimator(const std::string& name)
  : Estimator(name,"","scalar"), scale(1.), shift(0.) {
}

ScalarEstimator::ScalarEstimator(const std::string &name,
  const std::string &typeString, const std::string &unitName,
  double scale, double shift)
  : Estimator(name,typeString,unitName), scale(scale), shift(shift) {
}

void ScalarEstimator::averageOverClones(const MPIManager* mpi) {
  int rank=0, size=1;
#ifdef ENABLE_MPI
  if (mpi->isCloneMain()) {
    rank = mpi->getCloneComm().Get_rank();
    size = mpi->getCloneComm().Get_size();
  }
#endif
  double v=calcValue(),value=v;
  reset();
#ifdef ENABLE_MPI
  if (mpi->isCloneMain()) {
    mpi->getCloneComm().Reduce(&v,&value,1,MPI::DOUBLE,MPI::SUM,0);
  }
#endif
  if (rank==0) setValue(value/size);
}
