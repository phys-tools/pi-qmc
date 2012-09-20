#ifndef ASCIIPARTITIONEDSCALARREPORTWRITER_H_
#define ASCIIPARTITIONEDSCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
class ScalarEstimator;
class PartitionedScalarAccumulator;
#include <fstream>

class AsciiPartitionedScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, PartitionedScalarAccumulator> {
public:
    AsciiPartitionedScalarReportWriter(std::ofstream& file);
    virtual ~AsciiPartitionedScalarReportWriter();

    virtual void startReport(const ScalarEstimator* est,
            const PartitionedScalarAccumulator *acc);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator* est,
            const PartitionedScalarAccumulator *acc);

private:
    std::ofstream &file;
    int partitionCount;
};

#endif
