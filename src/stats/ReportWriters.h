#ifndef REPORTWRITERS_H_
#define REPORTWRITERS_H_

template <class T1, class T2> class ReportWriterInterface;
class ScalarEstimator;
class ArrayEstimator;
class AccRejEstimator;
class ScalarAccumulator;
class SimpleScalarAccumulator;
class PartitionedScalarAccumulator;

class ReportWriters {
public:
    ReportWriters(
            ReportWriterInterface<ScalarEstimator, SimpleScalarAccumulator>*,
            ReportWriterInterface<ScalarEstimator, PartitionedScalarAccumulator>*,
            ReportWriterInterface<ArrayEstimator, ScalarAccumulator>*,
            ReportWriterInterface<AccRejEstimator, ScalarAccumulator>*);
    virtual ~ReportWriters();

    void startScalarReport(ScalarEstimator *est, SimpleScalarAccumulator *acc);
    void reportScalarStep(ScalarEstimator *est, SimpleScalarAccumulator *acc);

    void startScalarReport(ScalarEstimator *est,
            PartitionedScalarAccumulator *acc);
    void reportScalarStep(ScalarEstimator *est,
            PartitionedScalarAccumulator *acc);

    void startAccRejReport(AccRejEstimator *est, ScalarAccumulator *acc);
    void reportAccRejStep(AccRejEstimator *est, ScalarAccumulator *acc);

    void startArrayReport(ArrayEstimator *est, ScalarAccumulator *acc);
    void reportArrayStep(ArrayEstimator *est, ScalarAccumulator *acc);

private:
    ReportWriterInterface<ScalarEstimator, SimpleScalarAccumulator>
        *scalarReportWriter;
    ReportWriterInterface<ScalarEstimator, PartitionedScalarAccumulator>
        *partitionedScalarReportWriter;
    ReportWriterInterface<ArrayEstimator, ScalarAccumulator>
        *arrayReportWriter;
    ReportWriterInterface<AccRejEstimator, ScalarAccumulator>
        *accRejReportWriter;
};

#endif
