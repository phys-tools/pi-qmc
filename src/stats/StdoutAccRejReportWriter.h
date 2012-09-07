#ifndef STDOUTACCREJREPORTWRITER_H_
#define STDOUTACCREJREPORTWRITER_H_

#include "ReportWriterInterface.h"
#include "AccRejEstimator.h"

class StdoutAccRejReportWriter: public ReportWriterInterface<AccRejEstimator> {
public:
    StdoutAccRejReportWriter();
    virtual ~StdoutAccRejReportWriter();

    virtual void reportStep(const AccRejEstimator &est);
};

#endif
