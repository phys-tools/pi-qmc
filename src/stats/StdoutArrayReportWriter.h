#ifndef STDOUTARRAYREPORTWRITER_H_
#define STDOUTARRAYREPORTWRITER_H_

#include "ReportWriterInterface.h"
#include "ArrayEstimator.h"

class StdoutArrayReportWriter: public ReportWriterInterface<ArrayEstimator> {
public:
    StdoutArrayReportWriter();
    virtual ~StdoutArrayReportWriter();

    virtual void reportStep(const ArrayEstimator &est);
};

#endif
