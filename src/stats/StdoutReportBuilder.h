#ifndef __StdoutReportBuilder_h_
#define __StdoutReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include <string>
#include <blitz/array.h>

class ReportWriters;
class StdoutScalarReportWriter;
class StdoutArrayReportWriter;
class StdoutAccRejReportWriter;

/** Class for reporting estimators to standard output.
 @author John Shumway */
class StdoutReportBuilder: public EstimatorReportBuilder {
public:
    StdoutReportBuilder();
    virtual ~StdoutReportBuilder();

    virtual void initializeReport(EstimatorManager*);
    virtual void collectAndWriteDataBlock(EstimatorManager*);

    virtual void reportScalarStep(const ScalarEstimator& est);
    virtual void startAccRejReport(const AccRejEstimator& est) {}
    virtual void reportAccRejStep(const AccRejEstimator& est);
    virtual void reportArrayBlockedStep(const ArrayEstimator& est);

private:
    int istep;
    int nstep;

    ReportWriters *reportWriters;
    StdoutScalarReportWriter *scalarWriter;
    StdoutArrayReportWriter *arrayWriter;
    StdoutAccRejReportWriter *accrejWriter;
    void writeBlockHeader();
};
#endif
