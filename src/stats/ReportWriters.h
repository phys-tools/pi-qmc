#ifndef REPORTWRITERS_H_
#define REPORTWRITERS_H_

class ScalarReportWriter;
class ArrayReportWriter;
class ArrayBlockedReportWriter;
class AccRejReportWriter;

class ReportWriters {
public:
    ReportWriters(ScalarReportWriter*, ArrayReportWriter*,
            ArrayBlockedReportWriter*, AccRejReportWriter*);
    virtual ~ReportWriters();

    ScalarReportWriter* getScalarReportWriter() const;
    ArrayReportWriter* getArrayReportWriter() const;
    ArrayBlockedReportWriter* getArrayBlockedReportWriter() const;
    AccRejReportWriter* getAccRejReportWriter() const;

private:
    ScalarReportWriter *scalarReportWriter;
    ArrayReportWriter *arrayReportWriter;
    ArrayBlockedReportWriter *arrayBlockedReportWriter;
    AccRejReportWriter *accRejReportWriter;
};

#endif
