#ifndef NULLPARTITIONEDSCALARREPORTWRITER_H_
#define NULLPARTITIONEDSCALARREPORTWRITER_H_

#include "ReportWriterInterface.h"
class ScalarEstimator;
class PartitionedScalarAccumulator;

class NullPartitionedScalrReportWriter
:   public ReportWriterInterface<ScalarEstimator, PartitionedScalarAccumulator> {
public:
    NullPartitionedScalrReportWriter() {}
    virtual ~NullPartitionedScalrReportWriter() {}

    virtual void startReport(const ScalarEstimator* estimator,
            const PartitionedScalarAccumulator* accumulator) {}
    virtual void startBlock(int istep) {}
    virtual void reportStep(const ScalarEstimator* estimator,
            const PartitionedScalarAccumulator* accumulator) {}
};

#endif
