#ifndef __ThermoEnergyEstimator_h_
#define __ThermoEnergyEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
class Paths;
class Action;
class DoubleAction;
class SimulationInfo;
class MPIManager;
/** Thermodynamic energy estimator. 
 * The thermodynamic energy estimator is obtained by differentiating
 *  the partition function with repect to @f$\beta@f$,
 *  @f[ E_T = -\frac{1}{Z}\frac{dZ}{d\beta} 
 *          = \left\langle \frac{dU_i}{d\tau} \right\rangle, @f]
 *  which is an average of the time derivative of the action.
 *  Note that Ceperley separates out the kinetic action in his
 *  review article (http://link.aps.org/abstract/RMP/v67/p279),
 *  we don't make this distinction. The time derivatives are obtained
 *  by calling Action::getBeadAction.
 *  @version $Revision$
 *  @author John Shumway  */
class ThermoEnergyEstimator : public ScalarEstimator, public LinkSummable {
public:
  /// Constructor.
  ThermoEnergyEstimator(const SimulationInfo& simInfo, const Action*,
          const DoubleAction*, MPIManager *mpi,
          const std::string& unitName, double scale, double shift);
  /// Virtual destructor.
  virtual ~ThermoEnergyEstimator() {}
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
  /// ThermoEnergy energy.
  double energy;
  /// The accumulated energy.
  double etot;
  /// The normalization of the accumualted energy.
  double enorm;
  /// Action.
  const Action* action;
  /// DoubleAction.
  const DoubleAction* doubleAction;
  /// The MPI Manager.
  MPIManager *mpi; 
};

#endif
