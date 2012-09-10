#ifndef STDOUTACCREJREPORTWRITER_H_
#define STDOUTACCREJREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
#include "stats/AccRejEstimator.h"

class StdoutAccRejReportWriter: public ReportWriterInterface<AccRejEstimator> {
public:
    StdoutAccRejReportWriter();
    virtual ~StdoutAccRejReportWriter();

    virtual void reportStep(const AccRejEstimator &est);
};

#endif
