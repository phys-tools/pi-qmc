#ifndef SCALARACUMULATOR_H_
#define SCALARACUMULATOR_H_

class MPIManager;
class ReportWriters;
class ScalarEstimator;

class ScalarAccumulator {
public:
    ScalarAccumulator(MPIManager *mpi);
    virtual ~ScalarAccumulator();
    void addToValue(double addend);
    void clearValue();
    void storeValue(const int lnslice);
    void reset();
    double calcValue();

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
