#include "AsciiPartitionedScalarReportWriter.h"
#include "stats/ScalarEstimator.h"
#include "stats/PartitionedScalarAccumulator.h"

AsciiPartitionedScalarReportWriter::AsciiPartitionedScalarReportWriter(std::ofstream &file)
:   file(file) {
}

AsciiPartitionedScalarReportWriter::~AsciiPartitionedScalarReportWriter() {
}

void AsciiPartitionedScalarReportWriter::startReport(const ScalarEstimator *est,
        const PartitionedScalarAccumulator *acc) {
    partitionCount = acc->getPartitionCount();
    for (int partition = 0; partition < partitionCount; ++ partition) {
        file << "\"" << est->getName() << partition;
        const std::string& unitName(est->getUnitName());
        if (unitName != "") {
            file << " (" << unitName << ")";
        }
        file << "\" ";
    }
}

void AsciiPartitionedScalarReportWriter::startBlock(int istep) {
}

void AsciiPartitionedScalarReportWriter::reportStep(const ScalarEstimator *est,
        const PartitionedScalarAccumulator *acc) {
    for (int partition = 0; partition < partitionCount; ++ partition) {
        double value
            = (acc->getValue(partition) + est->getShift()) * est->getScale();
        file << value << " ";
    }
}


