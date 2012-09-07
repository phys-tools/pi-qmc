#ifndef REPORTWRITERS_H_
#define REPORTWRITERS_H_

template <class T> class ReportWriterInterface;
class ScalarEstimator;
class ArrayEstimator;
class AccRejEstimator;

class ReportWriters {
public:
    ReportWriters(ReportWriterInterface<ScalarEstimator>*,
            ReportWriterInterface<ArrayEstimator>*,
            ReportWriterInterface<AccRejEstimator>*);
    virtual ~ReportWriters();

    void startScalarReport(ScalarEstimator &est);
    void reportScalarStep(ScalarEstimator &est);

    void startAccRejReport(AccRejEstimator &est);
    void reportAccRejStep(AccRejEstimator &est);

    void startArrayReport(ArrayEstimator &est);
    void reportArrayStep(ArrayEstimator &est);

private:
    ReportWriterInterface<ScalarEstimator> *scalarReportWriter;
    ReportWriterInterface<ArrayEstimator> *arrayReportWriter;
    ReportWriterInterface<AccRejEstimator> *accRejReportWriter;
};

#endif
