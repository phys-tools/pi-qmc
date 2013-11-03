#ifndef ARRAYACUMULATOR_H_
#define ARRAYACUMULATOR_H_

class ReportWriters;
class ArrayEstimator;

class ArrayAccumulator {
public:
    ArrayAccumulator() {}
    virtual ~ArrayAccumulator() {}
    virtual void addToValue(double addend) = 0;
    virtual void clearValue() = 0;
    virtual void storeValue(const int lnslice) = 0;
    virtual void averageOverClones() = 0;
    virtual void calculateTotal() = 0;
    virtual void reset() = 0;

    virtual void startReport(ReportWriters* writers,
            ArrayEstimator* estimator) = 0;
    virtual void reportStep(ReportWriters* writers,
            ArrayEstimator* estimator) = 0;
};

#endif
