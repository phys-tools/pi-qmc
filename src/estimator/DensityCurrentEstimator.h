#ifndef __DensityCurrentEstimator_h_
#define __DensityCurrentEstimator_h_
#include "stats/BlitzArrayBlkdEst.h"
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include <fftw3.h>
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <vector>
class Distance;
class Paths;
class SimulationInfo;
class MPIManager;
/// Density-current estimator estimator for an inhomogeneous system.
///
/// @f[
/// \chi_{nj}(x,\tau) = -\langle n(x,\tau) j_x(0,0) \rangle_beta
/// @f]
class DensityCurrentEstimator : public BlitzArrayBlkdEst<NDIM+2>, public LinkSummable {
public:
  typedef blitz::Array<std::complex<double>,NDIM+1> CArrayN;
  typedef blitz::Array<std::complex<double>,2> CArray2;
  typedef blitz::Array<std::complex<double>,1> CArray1;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<float,3> FArray3;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::TinyVector<int,NDIM+2> IVecN;
  typedef std::vector<Distance*> DistArray;

  /// Constructor.
  DensityCurrentEstimator(const SimulationInfo& simInfo,
    const std::string& name, const Vec &min, 
    const Vec &max, const IVec &nbin, const IVecN &nbinN, const int &njbin,
    const DistArray &dist, const int nstride, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~DensityCurrentEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const blitz::TinyVector<double,NDIM>& start,
                          const blitz::TinyVector<double,NDIM>& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  ///
  const int nsliceEff, nfreq, nstride, ntot;
  const Vec min;
  const Vec deltaInv;
  const IVec nbin;
  CArray2 tempj, tempj0;
  const double tau, ax;
  const DistArray dist;
  const int npart;
  Array q;
  const int njbin;
  const double dxinv;
  MPIManager *mpi;
  CArrayN tempn;
  CArray2* tempn_;
  FArray3* value_;
  fftw_plan fwdn, fwdj, fwdj0;
};

#endif
