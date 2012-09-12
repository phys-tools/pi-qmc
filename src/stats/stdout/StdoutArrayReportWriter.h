#ifndef STDOUTARRAYREPORTWRITER_H_
#define STDOUTARRAYREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
#include "stats/ArrayEstimator.h"
#include "stats/ScalarAccumulator.h"

class StdoutArrayReportWriter
:   public ReportWriterInterface<ArrayEstimator, ScalarAccumulator> {
public:
    StdoutArrayReportWriter();
    virtual ~StdoutArrayReportWriter();

    virtual void startReport(const ArrayEstimator* estimator,
            const ScalarAccumulator* accumulator);
    virtual void startBlock(int istep);
    virtual void reportStep(const ArrayEstimator* estimator,
            const ScalarAccumulator* accumulator);

};

#endif
