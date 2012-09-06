#ifndef __SKOmegaEstimator_h_
#define __SKOmegaEstimator_h_
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
#include <fftw3.h>

/** Calculates dynamic structure factor for a translationally invariant system.
 *  @version $Revision: 393 $
 *  @author John Shumway  */
class SKOmegaEstimator : public LinkSummable, 
                         public BlitzArrayBlkdEst<NDIM+3> {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef blitz::TinyVector<int, NDIM+3> IVecN;
  typedef blitz::Array<int,NDIM> IArrayNDIM;
  typedef blitz::Array<int,1> IArray;
  typedef std::complex<double> Complex;

  /// Constructor.
  SKOmegaEstimator(const SimulationInfo& simInfo, const std::string& name,
      const IVec &nbin, const IVecN &nbinN, int nstride, MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~SKOmegaEstimator();
  
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
  const int nslice,nfreq,nstride,npart,nspec,ntot;
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  SuperCell cell;
  const double tau;
  blitz::Array<Complex,NDIM+2> temp;
  IArray speciesIndex;
  fftw_plan fwd;
  MPIManager *mpi;
  blitz::Array<Complex,3>* temp_;
  blitz::Array<float,4>* value_;
};
#endif
