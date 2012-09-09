#ifndef STDOUTSCALARREPORTWRITER_H_
#define STDOUTSCALARREPORTWRITER_H_

#include "ReportWriterInterface.h"
#include "ScalarEstimator.h"
#include <blitz/array.h>

class StdoutScalarReportWriter: public ReportWriterInterface<ScalarEstimator> {
public:
    StdoutScalarReportWriter(int stepCount);
    virtual ~StdoutScalarReportWriter();

    virtual void startReport(const ScalarEstimator &est);
    virtual void reportStep(const ScalarEstimator &est);
    void startBlock(int istep);
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
