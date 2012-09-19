#ifndef PARTITIONEDCALARACUMULATOR_H_
#define PARTITIONEDSCALARACUMULATOR_H_

#include "ScalarAccumulator.h"
class MPIManager;
class PartitionWeight;

class PartitionedScalarAccumulator : public ScalarAccumulator {
public:
    PartitionedScalarAccumulator(MPIManager *mpi,
            PartitionWeight *modelState);
    virtual ~PartitionedScalarAccumulator();
    virtual void addToValue(double addend);
    virtual void clearValue();
    virtual void storeValue(const int lnslice);
    virtual void reset();
    double getValue(int partition) const;
    int getPartitionCount() const;

    virtual void startReport(ReportWriters* writers,
            ScalarEstimator* estimator);
    virtual void reportStep(ReportWriters* writers,
            ScalarEstimator* estimator);
private:
    double* value;
    double* sum;
    double* norm;
    const int partitionCount;
    MPIManager *mpi;
    PartitionWeight *weight;
};

#endif
