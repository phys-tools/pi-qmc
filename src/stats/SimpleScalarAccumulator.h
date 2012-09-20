#ifndef SIMPLESCALARACUMULATOR_H_
#define SIMPLESCALARACUMULATOR_H_

#include "ScalarAccumulator.h"
class MPIManager;

class SimpleScalarAccumulator : public ScalarAccumulator {
public:
    SimpleScalarAccumulator(MPIManager *mpi);
    virtual ~SimpleScalarAccumulator();
    virtual void addToValue(double addend);
    virtual void clearValue();
    virtual void storeValue(const int lnslice);
    virtual void reset();
    virtual void calculateTotal();
    double getValue() const;

    virtual void startReport(ReportWriters* writers,
            ScalarEstimator* estimator);
    virtual void reportStep(ReportWriters* writers,
            ScalarEstimator* estimator);
private:
    double value;
    double sum;
    double norm;
    MPIManager *mpi;
};

#endif
