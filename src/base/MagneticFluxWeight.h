#ifndef MAGNETICFLUXWEIGHT_H_
#define MAGNETICFLUXWEIGHT_H_

#include "stats/PartitionWeight.h"
class MagneticFluxCalculator;

class MagneticFluxWeight: public PartitionWeight {
public:
    MagneticFluxWeight(double bfield, int partitionCount,
            MagneticFluxCalculator*);
    virtual ~MagneticFluxWeight();
    virtual int getPartitionCount() const;
    virtual void evaluate(Paths* paths);
    virtual double getValue(int i) const;
private:
    double deltaB;
    double deltaPhi;
    int partitionCount;
    MagneticFluxCalculator* fluxCalculator;
};

#endif
