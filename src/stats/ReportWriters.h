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

    ReportWriterInterface<ScalarEstimator>* getScalarReportWriter() const;
    ReportWriterInterface<ArrayEstimator>* getArrayBlockedReportWriter() const;
    ReportWriterInterface<AccRejEstimator>* getAccRejReportWriter() const;

private:
    ReportWriterInterface<ScalarEstimator> *scalarReportWriter;
    ReportWriterInterface<ArrayEstimator> *arrayReportWriter;
    ReportWriterInterface<AccRejEstimator> *accRejReportWriter;
};

#endif
