#ifndef REPORTWRITERS_H_
#define REPORTWRITERS_H_

class ScalarReportWriter;
class ArrayReportWriter;
class AccRejReportWriter;

class ReportWriters {
public:
    ReportWriters(ScalarReportWriter*, ArrayReportWriter*,
            AccRejReportWriter*);
    virtual ~ReportWriters();

    ScalarReportWriter* getScalarReportWriter() const;
    ArrayReportWriter* getArrayBlockedReportWriter() const;
    AccRejReportWriter* getAccRejReportWriter() const;

private:
    ScalarReportWriter *scalarReportWriter;
    ArrayReportWriter *arrayReportWriter;
    AccRejReportWriter *accRejReportWriter;
};

#endif
