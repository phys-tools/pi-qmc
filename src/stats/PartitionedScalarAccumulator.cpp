#include "PartitionedScalarAccumulator.h"
#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MPIManager.h"
#include "ReportWriters.h"
#include "PartitionWeight.h"

PartitionedScalarAccumulator::PartitionedScalarAccumulator(MPIManager *mpi,
        PartitionWeight *weight)
:   partitionCount(weight->getPartitionCount()),
    value(new double[partitionCount]),
    sum(new double[partitionCount]),
    norm(new double[partitionCount]),
    mpi(mpi),
    weight(weight) {
}

PartitionedScalarAccumulator::~PartitionedScalarAccumulator() {
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
        double *buffer = new buffer[partitionCount];
        int ibuffer;
        mpi->getWorkerComm().Reduce(&value,&buffer,partitionCount,MPI::DOUBLE,MPI::SUM,0);
        mpi->getWorkerComm().Reduce(&lnslice,&ibuffer,1,MPI::INT,MPI::SUM,0);
        for (int i = 0; i < partitionCount; ++i) {
            value[i] = buffer[i];
        }
        nslice=ibuffer;
    }
#endif
    for (int i = 0; i < partitionCount; ++i) {
        value[i] /= nslice;
        sum[i] += value[i] * weight->getValue(i);
        norm[i] += 1.0;
    }
}

void PartitionedScalarAccumulator::reset() {
    for (int i = 0; i < partitionCount; ++i) {
        sum[i] = norm[i] = 0.0;
    }
}

void PartitionedScalarAccumulator::startReport(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->startScalarReport(estimator, this);
}

int PartitionedScalarAccumulator::getPartitionCount() const {
    return partitionCount;
}

double PartitionedScalarAccumulator::getValue(int partition) const {
    return sum[partition] / norm[partition];
}

void PartitionedScalarAccumulator::reportStep(ReportWriters* writers,
        ScalarEstimator* estimator) {
    writers->reportScalarStep(estimator, this);
}
