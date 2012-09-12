#ifndef STDOUTSCALARREPORTWRITER_H_
#define STDOUTSCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
#include "stats/ScalarEstimator.h"
#include "stats/ScalarAccumulator.h"
#include <blitz/array.h>

class StdoutScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, ScalarAccumulator> {
public:
    StdoutScalarReportWriter(int stepCount);
    virtual ~StdoutScalarReportWriter();

    virtual void startReport(const ScalarEstimator &est);
    virtual void reportStep(const ScalarEstimator &est);
    virtual void startBlock(int istep);
private:
    typedef blitz::Array<double, 1> Array;
    int nstep;
    int istep;
    int iscalar;
    Array sum;
    Array sum2;
    Array norm;
};

#endif
