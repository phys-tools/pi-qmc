#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "DensCountEstimator.h"
#include "action/Action.h"
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include "util/Distance.h"
#include <cstdlib>
#include <blitz/array.h>

DensCountEstimator::DensCountEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Species *spec,
    const Vec &min, const Vec &max, const IVec &nbin, const IVecN &nbinN,
    const DistArray &dist, MPIManager *mpi) 
  : BlitzArrayBlkdEst<NDIM+1>(name,"dynamic-array/density-count",nbinN,false),
    min(min), deltaInv(nbin/(max-min)), nbin(nbin), nbinN(nbinN), dist(dist),
    cell(*simInfo.getSuperCell()),temp(nbinN),count(nbin), 
    ifirst(spec->ifirst), npart(spec->count), mpi(mpi) {
  scale=new VecN(1.);
  origin=new VecN(0.);
  for (int i=0; i<NDIM; ++i) {
    (*scale)[i] = (max[i]-min[i])/nbin[i];
    (*origin)[i] = min[i];
  }
  value=0.;
}


DensCountEstimator::~DensCountEstimator() {
  for (int i=0; i<NDIM; ++i) delete dist[i];
  delete scale;
  delete origin;
}

void DensCountEstimator::initCalc(const int nslice,
    const int firstSlice) {
  temp=0;
  count=0;
}


void DensCountEstimator::handleLink(const Vec& start, const Vec& end,
    const int ipart, const int islice, const Paths &paths) {
  if (ipart>=ifirst && ipart<ifirst+npart) {
    Vec r=start;
    IVec ibin=0;
    for (int i=0; i<NDIM; ++i) {
      double d=(*dist[i])(r);
      ibin[i]=int(floor((d-min[i])*deltaInv[i]));
      if (ibin[i]<0 || ibin[i]>=nbin[i]) break;
      if (i==NDIM-1) ++count(ibin);
    }
  }
  if (ipart==ifirst+npart-1) {
    int ntot=count.size();
    int maxc=nbinN[NDIM];
    int *countPtr=count.data();
    float *tempPtr=temp.data();
    for (int i=0; i<ntot; ++i) {
      int c = *(countPtr+i);
      if (c<maxc) ++(*(tempPtr+maxc*i+c));
    }
    count=0;
  }
}


void DensCountEstimator::endCalc(const int nslice) {
  temp/=nslice;
  // First move all data to 1st worker. 
  int workerID=(mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
//    if (mpi) {
//      if (workerID==0) {
//        mpi->getWorkerComm().Reduce(&temp(0,0,0),&temp(0,0,0),
//                                    product(nbin),MPI::DOUBLE,MPI::SUM,0);
//      } else {
//        mpi->getWorkerComm().Reduce(MPI::IN_PLACE,&temp(0,0,0),
//                                    product(nbin),MPI::DOUBLE,MPI::SUM,0);
//      }
//    }
#endif
  ///Need code for multiple workers!
  if (workerID==0) {
    value+=temp;
    norm+=1;
  }
}

void DensCountEstimator::reset() {}

void DensCountEstimator::evaluate(const Paths& paths) {
  paths.sumOverLinks(*this);
}
