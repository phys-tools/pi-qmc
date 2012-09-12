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

    virtual void reportStep(const ArrayEstimator &est);
};

#endif
