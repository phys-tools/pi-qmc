#ifndef __DensityEstimator_h_
#define __DensityEstimator_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/BlitzArrayBlkdEst.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
class Distance;
/** Calculates single particle densities for many different geometries.
 *  Implements several options for studying fluctations.
 *  @version $Revision$
 *  @author John Shumway  */
class DensityEstimator : public LinkSummable, public BlitzArrayBlkdEst<NDIM> {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef std::vector<Distance*> DistArray;

  /// Constructor.
  DensityEstimator(const SimulationInfo& simInfo, const std::string& name,
      const Species *s, const Vec &min, const Vec &max,
      const IVec &nbin, const DistArray &dist, MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~DensityEstimator();
  
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);

  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
      const int ipart, const int islice, const Paths &paths);
  
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  
  /// Clear value of estimator.
  virtual void reset();

  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths);

protected:
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  const DistArray dist;
  int ifirst, npart;
  ArrayN temp;
private:
  SuperCell cell;
#ifdef ENABLE_MPI
  ArrayN mpiBuffer;
#endif
  MPIManager *mpi;
};
#endif
