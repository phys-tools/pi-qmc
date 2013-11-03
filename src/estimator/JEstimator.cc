#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "JEstimator.h"
#include "base/SimulationInfo.h"
#include "base/Paths.h"
#include "stats/MPIManager.h"
#include "util/Permutation.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <fstream>

JEstimator::JEstimator(const SimulationInfo& simInfo, int nbfield,
    double bmax, MPIManager *mpi)
  : BlitzArrayBlkdEst<2>("exchange", "array/singlet-triplet",
                         IVecN(2,nbfield), true), 
    npart(simInfo.getNPart()),
    mpi(mpi), bstep(bmax/(nbfield-1.)), nbfield(nbfield) {
  value = 0.;
  norm = 0;
  std::cout << "bstep,bmax = " << bstep << ", " << bmax << std::endl;
}

JEstimator::~JEstimator() {
}

void JEstimator::evaluate(const Paths &paths) {
  paths.sumOverLinks(*this);
  const Permutation &perm(paths.getPermutation());
  if (perm[0]==0) {
    value(0,0)+=1.;
    value(1,0)+=1.;
    for (int n=1; n<nbfield; ++n) {
      double phase=cos(area*n*bstep);
      value(0,n) += phase;
      value(1,n) += phase;
    }
  } else {
    value(0,0)-=1.;
    value(1,0)+=1.;
    for (int n=1; n<nbfield; ++n) {
      double phase=cos(area*n*bstep);
      value(0,n) -= phase;
      value(1,n) += phase;
    }
  }
  norm += 1.;
}

void JEstimator::initCalc(const int nslice, const int firstSlice) {
  area=0.;
};

void JEstimator::handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&) {
  area += 0.5*(start[0]*(end[1]-start[1])-start[1]*(end[0]-start[0]));
}

void JEstimator::endCalc(int nslice) {
}

