#include "SimpleScalarAccumulator.h"
#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MPIManager.h"
#include "ReportWriters.h"

SimpleScalarAccumulator::SimpleScalarAccumulator(MPIManager *mpi)
:   mpi(mpi) {
    reset();
}

SimpleScalarAccumulator::~SimpleScalarAccumulator() {
}

void SimpleScalarAccumulator::clearValue() {
    value = 0.0;
}

void SimpleScalarAccumulator::addToValue(double addend) {
    value += addend;
}

void SimpleScalarAccumulator::storeValue(const int lnslice) {
    int nslice = lnslice;
#ifdef ENABLE_MPI
    if (mpi) {
        double buffer; int ibuffer;
        mpi->getWorkerComm().Reduce(&value,&buffer,1,MPI::DOUBLE,MPI::SUM,0);
        mpi->getWorkerComm().Reduce(&lnslice,&ibuffer,1,MPI::INT,MPI::SUM,0);
        value=buffer; nslice=ibuffer;
    }
#endif
    value /= nslice;
    sum += value;
    norm += 1.0;
}

void SimpleScalarAccumulator::calculateTotal() {
    value = sum / norm;
}

void SimpleScalarAccumulator::reset() {
    sum = norm = 0.0;
}

double SimpleScalarAccumulator::getValue() const {
    return value;
}

void SimpleScalarAccumulator::startReport(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->startScalarReport(estimator, this);
}

void SimpleScalarAccumulator::averageOverClones() {
#ifdef ENABLE_MPI
    double v = value;
    mpi->getCloneComm().Reduce(&v,&value,1,MPI::DOUBLE,MPI::SUM,0);
    if (mpi->isCloneMain()) {
        value /= mpi->getNClone();
    }
#endif
}

void SimpleScalarAccumulator::reportStep(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->reportScalarStep(estimator, this);
}

