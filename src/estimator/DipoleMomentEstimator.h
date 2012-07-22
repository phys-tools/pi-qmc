#ifndef __DipoleMomentEstimator_h_
#define __DipoleMomentEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>
class Paths;
class SimulationInfo;
/** dipole moment estimator.
 *  @author Matthew Harowitz  */
class DipoleMomentEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor.
  DipoleMomentEstimator(const SimulationInfo& simInfo, const int index=2);
  /// Virtual destructor.
  virtual ~DipoleMomentEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of coulomb energy estimator.
  virtual double calcValue() {return value=dtot/dnorm;}
  /// Clear value of dipole energy estimator.
  virtual void reset() {dtot=dnorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// Dipole Moment.
  double dipole;
  /// The accumulated polarization.
  double dtot;
  /// The normalization of the accumualted polarization.
  double dnorm;
  /// The charges of the particles.
  Array q;
  /// The component of the dipole moment to measure (0->x, 1->y, 2->z)
  int index;
};

#endif
