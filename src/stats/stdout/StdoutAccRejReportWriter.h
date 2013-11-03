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

    virtual void startReport(const AccRejEstimator* estimator,
            const ScalarAccumulator* accumulator);
    virtual void startBlock(int istep);
    virtual void reportStep(const AccRejEstimator* estimator,
            const ScalarAccumulator* accumulator);

};

#endif
