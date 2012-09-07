#include "StdoutScalarReportWriter.h"

StdoutScalarReportWriter::StdoutScalarReportWriter(int stepCount,
        int estimatorCount) {
    nstep = stepCount;
    istep = 0;
    sum.resize(estimatorCount);
    sum2.resize(estimatorCount);
    norm.resize(estimatorCount);
    sum = 0;
    sum2 = 0;
    norm = 0;
}

StdoutScalarReportWriter::~StdoutScalarReportWriter() {
}

void StdoutScalarReportWriter::reportStep(const ScalarEstimator& est) {
    double value = est.getValue();
    sum(iscalar) += value;
    sum2(iscalar) += value * value;
    norm(iscalar) += 1;
    std::cout << est.getName();
    if (est.getUnitName() != "")
        std::cout << " (" << est.getUnitName() << ")";
    std::cout << ": " << value << ", " << "Av="
            << sum(iscalar) / (norm(iscalar)) << " +-"
            << sqrt(
                    sum2(iscalar)
                            - sum(iscalar) * sum(iscalar) / (norm(iscalar)))
                    / (norm(iscalar) - 1) << std::endl;
    ++iscalar;
}

void StdoutScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    iscalar = 0;
}
