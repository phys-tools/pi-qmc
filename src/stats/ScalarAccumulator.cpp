#include "ScalarAccumulator.h"
#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MPIManager.h"

ScalarAccumulator::ScalarAccumulator(MPIManager *mpi)
:   mpi(mpi) {
}

ScalarAccumulator::~ScalarAccumulator() {
}

void ScalarAccumulator::clearValue() {
    value = 0.0;
}

void ScalarAccumulator::addToValue(double addend) {
    value += addend;
}

void ScalarAccumulator::storeValue(const int lnslice) {
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

void ScalarAccumulator::reset() {
    sum = norm = 0.0;
}

double ScalarAccumulator::calcValue() {
    return sum / norm;
}


