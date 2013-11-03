#ifndef FERMIONWEIGHT_H_
#define FERMIONWEIGHT_H_

#include "stats/PartitionWeight.h"

class FermionWeight: public PartitionWeight {
    int sign;
public:
    FermionWeight();
    virtual ~FermionWeight();
    virtual int getPartitionCount() const;
    virtual void evaluate(Paths* paths);
    virtual double getValue(int i) const;
private:
};

#endif
