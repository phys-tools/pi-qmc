#ifndef STDOUTACCREJREPORTWRITER_H_
#define STDOUTACCREJREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
#include "stats/AccRejEstimator.h"
#include "stats/ScalarAccumulator.h"

class StdoutAccRejReportWriter
:   public ReportWriterInterface<AccRejEstimator, ScalarAccumulator> {
public:
    StdoutAccRejReportWriter();
    virtual ~StdoutAccRejReportWriter();

    virtual void reportStep(const AccRejEstimator &est);
};

#endif
