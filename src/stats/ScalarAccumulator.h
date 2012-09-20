#ifndef SCALARACUMULATOR_H_
#define SCALARACUMULATOR_H_

class ReportWriters;
class ScalarEstimator;

class ScalarAccumulator {
public:
    ScalarAccumulator() {}
    virtual ~ScalarAccumulator() {}
    virtual void addToValue(double addend) = 0;
    virtual void clearValue() = 0;
    virtual void storeValue(const int lnslice) = 0;
    virtual void calculateTotal() = 0;
    virtual void reset() = 0;

    virtual void startReport(ReportWriters* writers,
            ScalarEstimator* estimator) = 0;
    virtual void reportStep(ReportWriters* writers,
            ScalarEstimator* estimator) = 0;
};

#endif
