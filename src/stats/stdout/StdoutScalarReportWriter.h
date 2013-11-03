#ifndef STDOUTSCALARREPORTWRITER_H_
#define STDOUTSCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
class ScalarEstimator;
class SimpleScalarAccumulator;
#include <blitz/array.h>

class StdoutScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, SimpleScalarAccumulator> {
public:
    StdoutScalarReportWriter(int stepCount);
    virtual ~StdoutScalarReportWriter();

    virtual void startReport(const ScalarEstimator *est,
            const SimpleScalarAccumulator *acc);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator *est,
            const SimpleScalarAccumulator *acc);
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
