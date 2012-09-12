#include "StdoutArrayReportWriter.h"

StdoutArrayReportWriter::StdoutArrayReportWriter() {
}

StdoutArrayReportWriter::~StdoutArrayReportWriter() {
}

void StdoutArrayReportWriter::startReport(const ArrayEstimator* estimator,
        const ScalarAccumulator* accumulator) {
}

void StdoutArrayReportWriter::startBlock(int istep) {
}

void StdoutArrayReportWriter::reportStep(const ArrayEstimator *est,
        const ScalarAccumulator *acc) {
    std::cout << "(measured " << est->getName() << ")" << std::endl;
}
