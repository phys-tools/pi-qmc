#ifndef __DensDensEstimator_h_
#define __DensDensEstimator_h_
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
class Distance;
/** Calculates single particle densities for many different geometries.
 *  Implements several options for studying fluctations.
 *  @version $Revision$
 *  @author John Shumway  */
class DensDensEstimator : public LinkSummable, 
                          public BlitzArrayBlkdEst<2*NDIM+1> {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int, NDIM> IVec;
  typedef std::vector<Distance*> DistArray;
  typedef blitz::Array<int,NDIM> IArrayNDIM;
  typedef std::complex<double> Complex;

  /// Constructor.
  DensDensEstimator(const SimulationInfo& simInfo, const std::string& name,
      const Species *s1, const Species *s2,
      const Vec &min, const Vec &max, const IVec &nbin,
      const IVecN &nbinN, const DistArray &dist, int nstride,
      MPIManager *mpi); 

  /// Virtual destructor.
  virtual ~DensDensEstimator();
  
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
  const int nslice,nfreq,nstride,ntot;
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  const IVecN nbinN;
  const DistArray dist;
  SuperCell cell;
  const double tau;
  blitz::Array<Complex,NDIM+1> temp1;
  blitz::Array<Complex,NDIM+1> temp2;
  int ifirst1, npart1, ifirst2, npart2;
  fftw_plan fwd1, fwd2;
  MPIManager *mpi;
  blitz::Array<Complex,2>* temp1_;
  blitz::Array<Complex,2>* temp2_;
  blitz::Array<float,3>* value_;
};
#endif
