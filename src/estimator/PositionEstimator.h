#ifndef __PositionEstimator_h_
#define __PositionEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "stats/ScalarEstimator.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <iostream>
class Paths;
class SimulationInfo;
/** position estimator.
 *  @author Matthew Harowitz  */
class PositionEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<std::string,1> SArray;
  /// Constructor.
  PositionEstimator(const SimulationInfo& simInfo, const Species& species, const int index=2);
  /// Virtual destructor.
  virtual ~PositionEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of position estimator.
  virtual double calcValue() {return value=ptot/pnorm;}
  /// Clear value of position estimator.
  virtual void reset() {ptot=pnorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// Position .
  double position;
  /// The accumulated polarization.
  double ptot;
  /// The normalization of the accumualted polarization.
  double pnorm;
  /// The component of the position to measure (0->x, 1->y, 2->z)
  int index;
  /// The names of the particles.
  SArray names;
  /// The species to measure
  std::string spec;
};

#endif
