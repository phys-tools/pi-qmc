#ifndef __WindingEstimator_h_
#define __WindingEstimator_h_
#include "base/LinkSummable.h"
#include "base/Species.h"
#include "stats/BlitzArrayBlkdEst.h"
#include <blitz/array.h>
class Paths;
class SuperCell;
class SimulationInfo;
class MPIManager;
/** Estimator for winding around periodic boundary conditions.
 *  @version $Revision$
 *  @author John Shumway  */
class WindingEstimator : public BlitzArrayBlkdEst<NDIM>, public LinkSummable {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  WindingEstimator(const SimulationInfo& simInfo, int nmax,
    const std::string &name, bool isChargeCoupled, const Species *species, MPIManager *mpi);
  /// Virtual destructor.
  virtual ~WindingEstimator();
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
  int nmax;
  int npart;
  int ifirst;
  Vec winding;
  SuperCell &cell;
  IArray charge;
  bool isChargeCoupled;
  MPIManager *mpi;
};

#endif
