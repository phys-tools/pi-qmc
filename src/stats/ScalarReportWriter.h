#ifndef SCALARREPORTWRITER_H_
#define SCALARREPORTWRITER_H_

class ScalarEstimator;

class ScalarReportWriter {
public:
    ScalarReportWriter();
    virtual ~ScalarReportWriter();

    virtual void startScalarReport(const ScalarEstimator& est);
    virtual void reportScalarStep(const ScalarEstimator& est);
};

#endif
