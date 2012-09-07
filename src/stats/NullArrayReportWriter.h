#ifndef NULLARRAYREPORTWRITER_H_
#define NULLARRAYREPORTWRITER_H_

#include "ReportWriterInterface.h"
class ArrayEstimator;

class NullArrayReportWriter: public ReportWriterInterface<ArrayEstimator> {
public:
    NullArrayReportWriter() {}
    virtual ~NullArrayReportWriter(){}
};

#endif
