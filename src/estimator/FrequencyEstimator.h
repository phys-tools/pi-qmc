#ifndef __FrequencyEstimator_h_
#define __FrequencyEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/BlitzArrayBlkdEst.h"
#include <cstdlib>
#include <blitz/array.h>
#include <string.h>
#include <iostream>
#include <fftw3.h>
class Paths;
class Action;
class SimulationInfo;
class Species;
class MPIManager;
/** Frequency estimator.
 *  @version $Revision $ 
 *  @author Daejin Shin, John Shumway  */
class FrequencyEstimator : public BlitzArrayBlkdEst<1>,
                           public LinkSummable {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<std::complex<double>,1> CArray1;
  typedef blitz::Array<std::string,1> SArray;
  /// Constructor.
  FrequencyEstimator(const SimulationInfo& simInfo,
                     const Species& species1, const Species& species2,
                     int nfreq, int nstride, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~FrequencyEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int nfirst);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Clear value of bond-length estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// parameters.
  const int npart, nslice, nfreq, nstride;
  const double tau;
  /// Particle indices for the bond we are measuring.
  const int ipart, jpart;
  /// Temporary array for collecting correlation function. 
  CArray1 temp;
  /// FFT plan
  fftw_plan fwd;
  MPIManager *mpi;
#ifdef ENABLE_MPI
  CArray1 mpiBuffer;
#endif
};

#endif
