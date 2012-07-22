#ifndef __AngularMomentumEstimator_h_
#define __AngularMomentumEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
#include <string.h>
#include <math.h>
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
class Paths;
class PhaseModel;
class SimulationInfo;
/** Angular Momentum estimator.
 * Hard code for 2-dimensional systems.
 * @f[ L_z = \hbar (x \partial_y\phi_T - y \partial_x\phi_T) @f]
 *  @version $Revision$
 *  @author Daejin Shin  */
class AngularMomentumEstimator : public ScalarEstimator, public LinkSummable {
public:
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  AngularMomentumEstimator(const SimulationInfo& simInfo, PhaseModel*);
  /// Virtual destructor.
  virtual ~AngularMomentumEstimator();
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of coulomb energy estimator.
  virtual double calcValue() {return value=angMtot/angMnorm;}
  /// Clear value of coulomb energy estimator.
  virtual void reset() {angMtot=angMnorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// AngularMomentum energy.
  double angM;
  /// The accumulated energy.
  double angMtot;
  /// The normalization of the accumualted energy.
  double angMnorm;
  /// The phase model.
  PhaseModel* phaseModel;
  /// The number of particles.
  const int npart;
  /// The number of slices.
  const int nslice;
  /// Arrays for coordinates and phase gradient.
  VArray r1, r2, grad;
};


#endif
