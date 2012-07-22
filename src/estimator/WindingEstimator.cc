#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "WindingEstimator.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <fstream>

WindingEstimator::WindingEstimator(const SimulationInfo& simInfo,
  int nmax, const std::string &name, bool isChargeCoupled, const Species *species, MPIManager *mpi)
  : BlitzArrayBlkdEst<NDIM>(name,"histogram/winding", IVecN(2*nmax+1), true), 
    nmax(nmax), npart(simInfo.getNPart()), cell(*simInfo.getSuperCell()),
    charge(npart), isChargeCoupled(isChargeCoupled), ifirst(0), mpi(mpi) {
  value = 0.;
  norm = 0;
  charge = 1;
  if(species!=0) {
    ifirst = species->ifirst;
    npart = species->count;
  }

  if (isChargeCoupled) {
    for (int i=ifirst; i<npart; ++i) {
      charge(i)=int(simInfo.getPartSpecies(i).charge);
    }
  }
}

WindingEstimator::~WindingEstimator() {
}

void WindingEstimator::evaluate(const Paths &paths) {
  paths.sumOverLinks(*this);
}

void WindingEstimator::initCalc(const int nslice, const int firstSlice) {
  winding = 0.;
};

void WindingEstimator::handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec delta = start-end;
    cell.pbc(delta);
    winding += delta*charge(ipart);
  }
}

void WindingEstimator::endCalc(int nslice) {
  #ifdef ENABLE_MPI
  if (mpi) {
    Vec buffer;
    mpi->getWorkerComm().Reduce(&winding,&buffer,NDIM,MPI::DOUBLE,MPI::SUM,0);
    winding = buffer;
  }
  #endif
  IVec iwind = rint(winding*cell.b)+nmax;
  for (int idim=0; idim<NDIM; ++idim) {
    if (iwind[idim]<0 || iwind[idim]>2*nmax+1) break;
    if (idim==NDIM-1) ++value(iwind);
  }
  ++norm;
}

