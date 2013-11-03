#ifndef __JEstimator_h_
#define __JEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/BlitzArrayBlkdEst.h"
#include <fftw3.h>
class Paths;
class SimulationInfo;
class MPIManager;
/** Estimator for singlet-triplet exchange splitting.
 *  @version $Revision$
 *  @bug Hard coded for two particles and no MPI.
 *  @author John Shumway  */
class JEstimator : public BlitzArrayBlkdEst<2>, public LinkSummable {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<bool,1> BArray;
  typedef blitz::Array<double,1> Array;
  /// Constructor.
  JEstimator(const SimulationInfo& simInfo, int nbfield,
             double bmax, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~JEstimator();
  /// Clear value of the estimator.
  virtual void reset() {}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths);
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          int ipart, int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(int nslice);
private:
  ///
  const int npart;
  MPIManager *mpi;
  const double bstep;
  const int nbfield;
  double area;
};

#endif
