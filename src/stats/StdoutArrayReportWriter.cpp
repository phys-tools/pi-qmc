#include "StdoutArrayReportWriter.h"

StdoutArrayReportWriter::StdoutArrayReportWriter() {
}

StdoutArrayReportWriter::~StdoutArrayReportWriter() {
}

void StdoutArrayReportWriter::reportStep(const ArrayEstimator& est) {
    std::cout << "(measured " << est.getName() << ")" << std::endl;
}
