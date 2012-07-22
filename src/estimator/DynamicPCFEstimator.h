#ifndef __DynamicPCFEstimator_h_
#define __DynamicPCFEstimator_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "action/Action.h"
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/ArrayBlockedEstimator.h"
#include "stats/BlitzArrayBlkdEst.h"
#include "stats/MPIManager.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
#include <fftw3.h>
class PairDistance;
/** Calculates single particle densities for many different geometries.
 *  Implements several options for studying fluctations.
 *  @version $Revision$
 *  @author John Shumway  */
class DynamicPCFEstimator : public LinkSummable, 
                          public BlitzArrayBlkdEst<3> {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef blitz::Array<int,NDIM> IArrayNDIM;
  typedef std::complex<double> Complex;

  /// Constructor.
  DynamicPCFEstimator(const SimulationInfo& simInfo, const std::string& name,
    const Species *s1, const Species *s2,
    double min, double max, int nbin, int nfreq, 
    int nstride, const PairDistance *dist, MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~DynamicPCFEstimator();
  
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

private:
  const int nslice,nfreq,nstride,nbin;
  const double min;
  const double deltaInv;
  const IVecN nbinN;
  const PairDistance *dist;
  SuperCell cell;
  const double tau;
  blitz::Array<Complex,2> temp;
  int ifirst,jfirst,nipart,njpart;
  fftw_plan fwd;
  MPIManager *mpi;
};
#endif
