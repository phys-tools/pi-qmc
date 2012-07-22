#ifndef __PermutationEstimator_h_
#define __PermutationEstimator_h_
#include "base/Paths.h"
#include "base/Species.h"
#include "stats/BlitzArrayBlkdEst.h"
#include <fftw3.h>
#include <cstdlib>
#include <blitz/array.h>
class Paths;
class SimulationInfo;
class MPIManager;
/** Permutation estimator for a homogeneous system.
 *  @version $Revision$
 *  @bug Hard coded for two particles and no MPI.
 *  @author John Shumway  */

class PermutationEstimator : public BlitzArrayBlkdEst<1> {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<bool,1> BArray;
  typedef blitz::Array<double,1> Array;
  /// Constructor.
  PermutationEstimator(const SimulationInfo& simInfo,  const std::string& name,
		       const Species &s1,MPIManager *mpi);
  /// Virtual destructor.
  virtual ~PermutationEstimator();
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths);
  virtual void averageOverClones(const MPIManager* mpi);

private:
  int ifirst, nipart;
  const int npart;
  bool *visited; 
  MPIManager *mpi;
};

#endif
