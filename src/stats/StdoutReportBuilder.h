#ifndef __StdoutReportBuilder_h_
#define __StdoutReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <blitz/array.h>

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
    virtual void reportArrayBlockedStep(const ArrayBlockedEstimator& est);
private:
    typedef blitz::Array<double, 1> Array;
    int nstep;
    int istep;
    int iscalar;
    /// The sum of the estimator.
    Array sum, sum2, norm;
};
#endif
