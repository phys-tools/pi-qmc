#ifndef __CoulombEnergyEstimator_h_
#define __CoulombEnergyEstimator_h_
#include "stats/ScalarEstimator.h"
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include <cstdlib>
#include <blitz/array.h>
class Paths;
class SimulationInfo;
class ScalarAccumulator;

/** Coulomb energy estimator. 
 *  @author John Shumway  */
class CoulombEnergyEstimator: public ScalarEstimator, public LinkSummable {
public:
    typedef blitz::Array<double, 1> Array;
    typedef blitz::TinyVector<double, NDIM> Vec;
    CoulombEnergyEstimator(const SimulationInfo& simInfo,
            const double epsilon, const std::string& unitName,
            double scale, double shift, ScalarAccumulator*);
    virtual ~CoulombEnergyEstimator();
    virtual void initCalc(const int nslice, const int firstSlice);
    virtual void handleLink(const Vec& start, const Vec& end, const int ipart,
            const int islice, const Paths&);
    virtual void endCalc(const int nslice);
    virtual double calcValue();
    virtual void reset();
    virtual void evaluate(const Paths& paths);
private:
    const double epsilon;
    Array q;
};

#endif
