#include "AsciiScalarReportWriter.h"
#include "stats/ScalarEstimator.h"
#include "stats/SimpleScalarAccumulator.h"

AsciiScalarReportWriter::AsciiScalarReportWriter(std::ofstream &file)
:   file(file) {
}

AsciiScalarReportWriter::~AsciiScalarReportWriter() {
}

void AsciiScalarReportWriter::startReport(const ScalarEstimator *est,
        const SimpleScalarAccumulator *acc) {
    file << "\"" << est->getName();
    const std::string& unitName(est->getUnitName());
    if (unitName != "") {
        file << " (" << unitName << ")";
    }
    file << "\" ";
}

void AsciiScalarReportWriter::startBlock(int istep) {
}

void AsciiScalarReportWriter::reportStep(const ScalarEstimator *est,
        const SimpleScalarAccumulator *acc) {
    double value = 0.0;
    if (acc == 0) {
        value = est->getValue();
    } else {
        value = (acc->getValue() + est->getShift()) * est->getScale();
    }
    file << value << " ";
}


