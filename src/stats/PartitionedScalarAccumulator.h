#ifndef PARTITIONEDCALARACUMULATOR_H_
#define PARTITIONEDSCALARACUMULATOR_H_

#include "ScalarAccumulator.h"
class MPIManager;
class ModelStateInterface;

class PartitionedScalarAccumulator : public ScalarAccumulator {
public:
    PartitionedScalarAccumulator(MPIManager *mpi,
            ModelStateInterface *modelState);
    virtual ~PartitionedScalarAccumulator();
    virtual void addToValue(double addend);
    virtual void clearValue();
    virtual void storeValue(const int lnslice);
    virtual void reset();
    virtual double calcValue();
    int getPartitionCount() const;

    virtual void startReport(ReportWriters* writers,
            ScalarEstimator* estimator);
    virtual void reportStep(ReportWriters* writers,
            ScalarEstimator* estimator);
private:
    double value;
    double sum;
    double norm;
    MPIManager *mpi;
    ModelStateInterface *modelState;
};

#endif
