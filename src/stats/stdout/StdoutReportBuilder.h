#ifndef __StdoutReportBuilder_h_
#define __StdoutReportBuilder_h_

#include "stats/EstimatorReportBuilder.h"
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
