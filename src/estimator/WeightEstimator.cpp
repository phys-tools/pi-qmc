#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "WeightEstimator.h"
#include "stats/ScalarAccumulator.h"

WeightEstimator::WeightEstimator(ScalarAccumulator *accumulator)
:   ScalarEstimator("weight", "scalar/weight", "", 1.0, 0.0) {
    this->accumulator = accumulator;
}

WeightEstimator::~WeightEstimator() {
}

double WeightEstimator::calcValue() {
    return 0.0;
}

void WeightEstimator::reset() {
    accumulator->reset();
}

void WeightEstimator::evaluate(const Paths& paths) {
    accumulator->clearValue();
    accumulator->addToValue(1.0);
    accumulator->storeValue(1);
}


