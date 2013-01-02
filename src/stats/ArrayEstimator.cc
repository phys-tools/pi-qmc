#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ArrayEstimator.h"
#include "ArrayAccumulator.h"

ArrayEstimator::ArrayEstimator(const std::string &name,
        const std::string &typeString, bool hasError)
    :   Estimator(name, typeString, ""),
        hasErrorFlag(hasError),
        accumulator(0) {
}

ArrayEstimator::~ArrayEstimator() {
    delete accumulator;
}

void ArrayEstimator::startReport(ReportWriters* writers) {
    writers->startArrayReport(this, 0);
}

void ArrayEstimator::reportStep(ReportWriters* writers) {
    writers->reportArrayStep(this, 0);
}

bool ArrayEstimator::hasError() const {
    return hasErrorFlag;
}



