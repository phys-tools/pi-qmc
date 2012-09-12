#ifndef NULLARRAYREPORTWRITER_H_
#define NULLARRAYREPORTWRITER_H_

#include "ReportWriterInterface.h"
class ArrayEstimator;
class ScalarAccumulator;

class NullArrayReportWriter
:   public ReportWriterInterface<ArrayEstimator, ScalarAccumulator> {
public:
    NullArrayReportWriter() {}
    virtual ~NullArrayReportWriter(){}
};

#endif
