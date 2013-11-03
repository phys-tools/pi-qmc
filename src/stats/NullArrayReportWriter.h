#ifndef NULLARRAYREPORTWRITER_H_
#define NULLARRAYREPORTWRITER_H_

#include "ReportWriterInterface.h"
class ArrayEstimator;
class ScalarAccumulator;

class NullArrayReportWriter
:   public ReportWriterInterface<ArrayEstimator, ScalarAccumulator> {
public:
    NullArrayReportWriter() {}
    virtual ~NullArrayReportWriter() {}

    virtual void startReport(const ArrayEstimator* estimator,
            const ScalarAccumulator* accumulator) {}
    virtual void startBlock(int istep) {}
    virtual void reportStep(const ArrayEstimator* estimator,
            const ScalarAccumulator* accumulator) {}
};

#endif
