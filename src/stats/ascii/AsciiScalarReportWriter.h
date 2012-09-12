#ifndef ASCIISCALARREPORTWRITER_H_
#define ASCIISCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
#include "stats/ScalarEstimator.h"
#include "stats/ScalarAccumulator.h"
#include <fstream>

class AsciiScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, ScalarAccumulator> {
public:
    AsciiScalarReportWriter(std::ofstream& file);
    virtual ~AsciiScalarReportWriter();

    virtual void startReport(const ScalarEstimator* est,
            const ScalarAccumulator *acc);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator* est,
            const ScalarAccumulator *acc);

private:
    std::ofstream &file;
};

#endif
