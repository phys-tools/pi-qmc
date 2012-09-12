#ifndef REPORTWRITERS_H_
#define REPORTWRITERS_H_

template <class T1, class T2> class ReportWriterInterface;
class ScalarEstimator;
class ArrayEstimator;
class AccRejEstimator;
class ScalarAccumulator;

class ReportWriters {
public:
    ReportWriters(ReportWriterInterface<ScalarEstimator, ScalarAccumulator>*,
            ReportWriterInterface<ArrayEstimator, ScalarAccumulator>*,
            ReportWriterInterface<AccRejEstimator, ScalarAccumulator>*);
    virtual ~ReportWriters();

    void startScalarReport(ScalarEstimator *est, ScalarAccumulator *acc);
    void reportScalarStep(ScalarEstimator *est, ScalarAccumulator *acc);

    void startAccRejReport(AccRejEstimator *est, ScalarAccumulator *acc);
    void reportAccRejStep(AccRejEstimator *est, ScalarAccumulator *acc);

    void startArrayReport(ArrayEstimator *est, ScalarAccumulator *acc);
    void reportArrayStep(ArrayEstimator *est, ScalarAccumulator *acc);

private:
    ReportWriterInterface<ScalarEstimator, ScalarAccumulator>
        *scalarReportWriter;
    ReportWriterInterface<ArrayEstimator, ScalarAccumulator>
        *arrayReportWriter;
    ReportWriterInterface<AccRejEstimator, ScalarAccumulator>
        *accRejReportWriter;
};

#endif
