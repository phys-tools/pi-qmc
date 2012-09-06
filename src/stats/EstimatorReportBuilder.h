#ifndef __EstimatorReportBuilder_h_
#define __EstimatorReportBuilder_h_

#include <string>
#include "ScalarReportWriter.h"
class EstimatorManager;
class AccRejEstimator;
class ArrayEstimator;

/** Base class for reporting estimators.
 @author John Shumway */

class EstimatorReportBuilder : public ScalarReportWriter {
public:
    EstimatorReportBuilder();
    virtual ~EstimatorReportBuilder();

    virtual void initializeReport(EstimatorManager*) = 0;
    virtual void collectAndWriteDataBlock(EstimatorManager*) = 0;

    virtual void startScalarReport(const ScalarEstimator& est);
    virtual void reportScalarStep(const ScalarEstimator& est);

    virtual void startArrayReport(const ArrayEstimator& est);
    virtual void reportArrayStep(const ArrayEstimator& est);

    virtual void startAccRejReport(const AccRejEstimator& est);
    virtual void reportAccRejStep(const AccRejEstimator& est);

    virtual void recordInputDocument(const std::string &docstring);
};
#endif
