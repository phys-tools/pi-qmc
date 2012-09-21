#include "PartitionedScalarAccumulator.h"
#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MPIManager.h"
#include "ReportWriters.h"
#include "PartitionWeight.h"
#include <iostream>

PartitionedScalarAccumulator::PartitionedScalarAccumulator(MPIManager *mpi,
        PartitionWeight *weight)
:   partitionCount(weight->getPartitionCount()),
    value(new double[partitionCount]),
    sum(new double[partitionCount]),
    norm(new double[partitionCount]),
    mpi(mpi),
    weight(weight) {
    reset();
}

PartitionedScalarAccumulator::~PartitionedScalarAccumulator() {
    delete [] value;
    delete [] sum;
    delete [] norm;
}

void PartitionedScalarAccumulator::clearValue() {
    for (int i = 0; i < partitionCount; ++i) {
        value[i] = 0.0;
    }
}

void PartitionedScalarAccumulator::addToValue(double addend) {
    for (int i = 0; i < partitionCount; ++i) {
        value[i] += addend;
    }
}

void PartitionedScalarAccumulator::storeValue(const int lnslice) {
    int nslice = lnslice;
#ifdef ENABLE_MPI
    if (mpi) {
        double *buffer = new double[partitionCount];
        int ibuffer;
        mpi->getWorkerComm().Reduce(value,buffer,partitionCount,MPI::DOUBLE,MPI::SUM,0);
        mpi->getWorkerComm().Reduce(&lnslice,&ibuffer,1,MPI::INT,MPI::SUM,0);
        for (int i = 0; i < partitionCount; ++i) {
            value[i] = buffer[i];
        }
        nslice=ibuffer;
        delete [] buffer;
    }
#endif
    for (int i = 0; i < partitionCount; ++i) {
        value[i] /= nslice;
        sum[i] += value[i] * weight->getValue(i);
        norm[i] += 1.0;
    }
}

void PartitionedScalarAccumulator::calculateTotal() {
    for (int i = 0; i < partitionCount; ++i) {
        value[i] = sum[i] / norm[i];
    }
}

void PartitionedScalarAccumulator::reset() {
    for (int i = 0; i < partitionCount; ++i) {
        sum[i] = norm[i] = 0.0;
    }
}

double PartitionedScalarAccumulator::getValue(int partition) const {
    return value[partition];
}

void PartitionedScalarAccumulator::startReport(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->startScalarReport(estimator, this);
}

int PartitionedScalarAccumulator::getPartitionCount() const {
    return partitionCount;
}

void PartitionedScalarAccumulator::averageOverClones() {
#ifdef ENABLE_MPI
    double *buffer = new double[partitionCount];
    for (int i = 0; i < partitionCount; ++i) {
        buffer[i] = value[i];
    }
    mpi->getCloneComm().Reduce(buffer,value,partitionCount,MPI::DOUBLE,MPI::SUM,0);
    delete [] buffer;
    if (mpi->isCloneMain()) {
        double oneOverSize = 1.0 / mpi->getNClone();
        for (int i = 0; i < partitionCount; ++i) {
            value[i] *= oneOverSize;
        }
    }
#endif
}

void PartitionedScalarAccumulator::reportStep(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->reportScalarStep(estimator, this);
}
