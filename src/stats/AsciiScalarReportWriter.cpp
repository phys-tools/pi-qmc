#include "AsciiScalarReportWriter.h"

AsciiScalarReportWriter::AsciiScalarReportWriter(std::ofstream &file)
:   file(file) {
}

AsciiScalarReportWriter::~AsciiScalarReportWriter() {
}

void AsciiScalarReportWriter::startReport(const ScalarEstimator& est) {
    file << "\"" << est.getName();
    const std::string& unitName(est.getUnitName());
    if (unitName != "") {
        file << " (" << unitName << ")";
    }
    file << "\" ";
}

void AsciiScalarReportWriter::reportStep(const ScalarEstimator& est) {
    file << est.getValue() << " ";
}


