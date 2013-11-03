#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "AccRejEstimator.h"
#include "MPIManager.h"

AccRejEstimator::AccRejEstimator(const std::string& name, const int nlevel)
  : Estimator(name,"acc-rej/multilevel",""),
    nlevel(nlevel), naccept(nlevel), ntrial(nlevel),
    sum(nlevel) {
  naccept=0; ntrial=0; sum=0;
}

void AccRejEstimator::averageOverClones(const MPIManager* mpi) {
#ifdef ENABLE_MPI
  MPI::COMM_WORLD.Reduce(naccept.data(),sum.data(),nlevel,MPI::LONG,MPI::SUM,0);
  naccept=0;
  if (mpi->isCloneMain()) naccept=sum;
  MPI::COMM_WORLD.Reduce(ntrial.data(),sum.data(),nlevel,MPI::LONG,MPI::SUM,0);
  ntrial=0;
  if (mpi->isCloneMain()) ntrial=sum;
#endif
}
