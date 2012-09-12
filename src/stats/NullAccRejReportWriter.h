#ifndef NULLACCREJREPORTWRITER_H_
#define NULLACCREJREPORTWRITER_H_

#include "ReportWriterInterface.h"
class AccRejEstimator;

class NullAccRejReportWriter
:   public ReportWriterInterface<AccRejEstimator, ScalarAccumulator> {
public:
    NullAccRejReportWriter() {}
    virtual ~NullAccRejReportWriter(){}
};

#endif
