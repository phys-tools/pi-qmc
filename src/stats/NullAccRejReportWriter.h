#ifndef NULLACCREJREPORTWRITER_H_
#define NULLACCREJREPORTWRITER_H_

#include "ReportWriterInterface.h"
#include "stats/AccRejEstimator.h"
#include "stats/ScalarAccumulator.h"

class NullAccRejReportWriter
:   public ReportWriterInterface<AccRejEstimator, ScalarAccumulator> {
public:
    NullAccRejReportWriter() {}
    virtual ~NullAccRejReportWriter() {}

    virtual void startReport(const AccRejEstimator* estimator,
            const ScalarAccumulator* accumulator) {}
    virtual void startBlock(int istep) {}
    virtual void reportStep(const AccRejEstimator* estimator,
            const ScalarAccumulator* accumulator) {}
};

#endif
