#ifndef __EMARateEstimator_h_
#define __EMARateEstimator_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/ScalarEstimator.h"
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>
class Paths;
class SimulationInfo;
class SuperCell;
class Species;
class CoulombLinkAction;

class EMARateEstimator : public ScalarEstimator, public LinkSummable {
public:
    typedef blitz::Array<double,1> Array;
    typedef blitz::TinyVector<double,NDIM> Vec;
    EMARateEstimator(const SimulationInfo& simInfo,
            const Species* species1, const Species* species2, double C);
    virtual ~EMARateEstimator();
    void includeCoulombContribution(double epsilon, int norder);

    virtual void initCalc(const int nslice, const int firstSlice);
    virtual void handleLink(const Vec& start, const Vec& end,
            const int ipart, const int islice, const Paths&);
    virtual void endCalc(const int nslice);
    virtual double calcValue() {return sum/norm;}
    virtual void reset() {sum=norm=0.;}
    virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
    const double dtau;
    Vec mass1;
    Vec mass2;
    const double C;
    const SuperCell& cell;
    const int lastSlice;
    double actionDifference;
    double sum, norm;
    bool hasCoulomb;
    CoulombLinkAction* coulomb;
    const Species* species1;
    const Species* species2;
    const int index1;
    const int index2;
    void evaluateElectronBeforeRecombination(
            const Vec& start, const Vec& end, const Paths & paths);
    void evaluateHoleAfterRecombination(
            const Vec& start, const Vec& end, const Paths & paths);

};

#endif
