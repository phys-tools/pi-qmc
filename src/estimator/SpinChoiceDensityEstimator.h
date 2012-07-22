#ifndef __SpinChoiceDensityEstimator_h_
#define __SpinChoiceDensityEstimator_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "DensityEstimator.h"
#include "base/LinkSummable.h"
#include "base/ModelState.h"
#include "base/Paths.h"
#include "base/SpinModelState.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/ArrayBlockedEstimator.h"
#include "stats/BlitzArrayBlkdEst.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include "util/Distance.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
/** Calculates single particle densities for many different geometries with
    SpinChoice action on.  */
class SpinChoiceDensityEstimator : public DensityEstimator {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef std::vector<Distance*> DistArray;

  /// Constructor.
  SpinChoiceDensityEstimator(const SimulationInfo& simInfo, const std::string& name,
      const Species *spec, const Vec &min, const Vec &max,
      const IVec &nbin, const DistArray &dist, int ispin, 
      ModelState& modelState, MPIManager *mpi) 
    : DensityEstimator(simInfo, name, spec, min, max, nbin, dist, mpi),
      ispin(ispin),
      spinState(dynamic_cast<SpinModelState&>(modelState).getSpinState()) {
  }

  /// Virtual destructor.
  virtual ~SpinChoiceDensityEstimator() {};
  
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
      const int ipart, const int islice, const Paths &paths) {
    if (spinState(ipart) != ispin) return;
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

private:
  int ispin;
  const IArray &spinState;
};
#endif
