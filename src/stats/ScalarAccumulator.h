#ifndef SCALARACUMULATOR_H_
#define SCALARACUMULATOR_H_

class MPIManager;

class ScalarAccumulator {
public:
    ScalarAccumulator(MPIManager *mpi);
    virtual ~ScalarAccumulator();
    void addToValue(double addend);
    void clearValue();
    void storeValue(const int lnslice);
    void reset();
    double calcValue();
private:
    double value;
    double sum;
    double norm;
    MPIManager *mpi;
};

#endif
