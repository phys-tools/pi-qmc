#ifndef __ThermoEnergyEstimator_h_
#define __ThermoEnergyEstimator_h_
#include "base/LinkSummable.h"
#include "stats/ScalarEstimator.h"
class Paths;
class Action;
class DoubleAction;
class SimulationInfo;
class ScalarAccumulator;
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
 *  @author John Shumway  */
class ThermoEnergyEstimator: public ScalarEstimator, public LinkSummable {
public:
    ThermoEnergyEstimator(const SimulationInfo& simInfo, const Action*,
            const DoubleAction*, const std::string& unitName,
            double scale, double shift, ScalarAccumulator*);
    virtual ~ThermoEnergyEstimator();
    virtual void initCalc(const int nslice, const int firstSlice);
    virtual void handleLink(const blitz::TinyVector<double, NDIM>& start,
            const blitz::TinyVector<double, NDIM>& end, const int ipart,
            const int islice, const Paths&);
    virtual void endCalc(const int nslice);
    virtual double calcValue();
    virtual void reset();
    virtual void evaluate(const Paths& paths);
private:
    ScalarAccumulator *accumulator;
    const Action* action;
    const DoubleAction* doubleAction;
};

#endif
