#include "EstimatorIterator.h"

EstimatorIterator::EstimatorIterator(
        std::list<Estimator*>& list)
:   iterator(list.begin()),
    end(list.end()) {
}

EstimatorIterator::~EstimatorIterator() {
}

bool EstimatorIterator::step() {
    return ++iterator != end;
}

Estimator* EstimatorIterator::operator*() {
    return *iterator;
}



