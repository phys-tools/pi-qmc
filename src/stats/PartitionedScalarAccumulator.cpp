#include "PartitionedScalarAccumulator.h"
#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MPIManager.h"
#include "ReportWriters.h"

#include <iostream>

PartitionedScalarAccumulator::PartitionedScalarAccumulator(MPIManager *mpi,
        ModelStateInterface *modelState)
:   mpi(mpi),
    modelState(modelState) {
}

PartitionedScalarAccumulator::~PartitionedScalarAccumulator() {
}

void PartitionedScalarAccumulator::clearValue() {
    value = 0.0;
}

void PartitionedScalarAccumulator::addToValue(double addend) {
    value += addend;
}

void PartitionedScalarAccumulator::storeValue(const int lnslice) {
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

void PartitionedScalarAccumulator::reset() {
    sum = norm = 0.0;
}

double PartitionedScalarAccumulator::calcValue() {
    return sum / norm;
}

void PartitionedScalarAccumulator::startReport(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->startScalarReport(estimator, this);
}

void PartitionedScalarAccumulator::reportStep(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->reportScalarStep(estimator, this);
}



