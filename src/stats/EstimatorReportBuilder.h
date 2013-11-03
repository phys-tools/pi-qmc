#ifndef __EstimatorReportBuilder_h_
#define __EstimatorReportBuilder_h_

#include <string>
class EstimatorManager;

/** Base class for reporting estimators.
 @author John Shumway */
class EstimatorReportBuilder {
public:
    EstimatorReportBuilder();
    virtual ~EstimatorReportBuilder();

    virtual void initializeReport(EstimatorManager*) = 0;
    virtual void collectAndWriteDataBlock(EstimatorManager*) = 0;
    virtual void recordInputDocument(const std::string &docstring);
};
#endif
