#ifndef __BoxEstimator_h_
#define __BoxEstimator_h_
#include "base/LinkSummable.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include <iostream>
class Paths;
class SimulationInfo;
/** box estimator.
 *  @todo add 3-d box and sphere
 *  @todo add 2-d box and circle
 *  @author Matthew Harowitz  */
class BoxEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<std::string,1> SArray;
  /// Constructor.
  BoxEstimator(const SimulationInfo& simInfo, const Species& species,
               const double lower, const double upper, const int index=2);
  /// Virtual destructor.
  virtual ~BoxEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of estimator.
  virtual double calcValue() {return value=intot/innorm;}
  /// Clear value of estimator.
  virtual void reset() {intot=innorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// In box?.
  double in;
  /// The accumulated in.
  double intot;
  /// The normalization of in.
  int innorm;
  /// The lower bound
  double low;
  /// The upper bound
  double up;
  /// The component of the box to measure (0->x, 1->y, 2->z)
  int index;
  /// The names of the particles.
  SArray names;
  /// The species to measure
  std::string spec;
};

#endif
