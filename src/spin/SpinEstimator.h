#ifndef __SpinEstimator_h_
#define __SpinEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
class Paths;
class Action;
class DoubleAction;
class SimulationInfo;
#include <string>
/** Spin state estimator.
 * @f[   \langle \uparrow | \sigma_{z} |\uparrow \rangle 
 *	 = \int d^{4}s \chi^{*}_{\downarrow}(\vec s) V(s) 
 * 	 \chi(\vec s)_{\uparrow}(s) = 1  @f]
 *  where @f[ V(\vec s)=6(1+2\alpha)^3 \dfrac{s_{0}^{2}-s_{1}^{2}-s_{2}^{2}+s_{3}^{2}}{2}
 *                 \dfrac{exp{-\alpha s^2}}{s^2} @f]
 *  
 *   @version $Revision$
 **  @author John Shumway, Daejin Shin  */
class SpinEstimator : public ScalarEstimator, public LinkSummable {
public:
  /// Constructor.
  SpinEstimator(const SimulationInfo& simInfo, const int idim, const double gc);
  /// Virtual destructor.
  virtual ~SpinEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const blitz::TinyVector<double,NDIM>& start,
                          const blitz::TinyVector<double,NDIM>& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of kinetic energy estimator.
  virtual double calcValue() {return value=etot/enorm;}
  /// Clear value of the estimator.
  virtual void reset() {etot=enorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// The dimension.
  const int idim;
  /// The Gaussian coefficient.
  const double gc, invgc; 
  /// Spin energy.
  double energy;
  /// The accumulated energy.
  double etot;
  /// The normalization of the accumualted energy.
  double enorm;
  /// Action.
  const Action* action;
  /// DoubleAction.
  const DoubleAction* doubleAction;
  /// Spin names.
  static const std::string namex,namey,namez;
  /// The width of the spin gaussian (l2=1/(mass*omega)).
  const double l2;
};

#endif
