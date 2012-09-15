#ifndef ASCIISCALARREPORTWRITER_H_
#define ASCIISCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
class ScalarEstimator;
class SimpleScalarAccumulator;
#include <fstream>

class AsciiScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, SimpleScalarAccumulator> {
public:
    AsciiScalarReportWriter(std::ofstream& file);
    virtual ~AsciiScalarReportWriter();

    virtual void startReport(const ScalarEstimator* est,
            const SimpleScalarAccumulator *acc);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator* est,
            const SimpleScalarAccumulator *acc);

private:
    std::ofstream &file;
};

#endif
