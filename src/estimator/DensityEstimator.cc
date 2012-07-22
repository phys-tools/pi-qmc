#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "DensityEstimator.h"
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include "util/Distance.h"
#include <cstdlib>
#include <blitz/array.h>

DensityEstimator::DensityEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species *spec,
    const Vec &min, const Vec &max, const IVec &nbin,
    const DistArray &dist, MPIManager *mpi) 
  : BlitzArrayBlkdEst<NDIM>(name,"array/density",nbin,false),
    min(min), deltaInv(nbin/(max-min)), nbin(nbin), dist(dist),
    ifirst(spec->ifirst), npart(spec->count), temp(nbin),
    cell(*simInfo.getSuperCell()),
#ifdef ENABLE_MPI
    mpiBuffer(nbin),
#endif 
    mpi(mpi) {
  scale=new Vec((max-min)/nbin);
  origin=new Vec(min);
  value=0.;
}


DensityEstimator::~DensityEstimator() {
  for (int i=0; i<NDIM; ++i) delete dist[i];
  delete scale;
  delete origin;
}

void DensityEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0.;
}


void DensityEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec r=start;
    IVec ibin=0;
    for (int i=0; i<NDIM; ++i) {
      double d=(*dist[i])(r);
      ibin[i]=int(floor((d-min[i])*deltaInv[i]));
      if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
      if (i==NDIM-1) ++temp(ibin);
    }
  }
}


void DensityEstimator::endCalc(const int lnslice) {
  int nslice = lnslice;
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    int ibuffer;
    mpi->getWorkerComm().Reduce(temp.data(),mpiBuffer.data(),
                                product(nbin),MPI::FLOAT,MPI::SUM,0);
    mpi->getWorkerComm().Reduce(&lnslice,&ibuffer,1,MPI::INT,MPI::SUM,0);
    temp = mpiBuffer;
    nslice = ibuffer;
  }
#endif
  temp /= nslice;
  if (workerID==0) {
    BlitzArrayBlkdEst<NDIM>::value+=temp;
    norm+=1.;
  }
}

void DensityEstimator::reset() {}

void DensityEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
