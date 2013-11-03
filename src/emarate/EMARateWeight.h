#ifndef __EMARateWeight_h_
#define __EMARateWeight_h_
#include "base/LinkSummable.h"
#include "base/Paths.h"
#include "stats/PartitionWeight.h"
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>
class Paths;
class SimulationInfo;
class SuperCell;
class Species;
class CoulombLinkAction;

class EMARateWeight : public PartitionWeight {
public:
    typedef blitz::Array<double,1> Array;
    typedef blitz::TinyVector<double,NDIM> Vec;
    EMARateWeight(const SimulationInfo& simInfo,
            const Species* species1, const Species* species2, double C);
    virtual ~EMARateWeight();
    void includeCoulombContribution(double epsilon, int norder);

    virtual int getPartitionCount() const;
    virtual void evaluate(Paths* paths);
    virtual double getValue(int i) const;
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
    double diagonalProbability;
    double recombiningProbability;
    void evaluateElectronBeforeRecombination(
            const Vec& start, const Vec& end, const Paths & paths);
    void evaluateHoleAfterRecombination(
            const Vec& start, const Vec& end, const Paths & paths);

};

#endif
