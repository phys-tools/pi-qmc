#include "ReportWriters.h"

ReportWriters::ReportWriters(ScalarReportWriter* scalarReportWriter,
        ArrayReportWriter* arrayBlockedReportWriter,
        AccRejReportWriter* accRejReportWriter)
:   scalarReportWriter(scalarReportWriter),
    arrayReportWriter(arrayBlockedReportWriter),
    accRejReportWriter(accRejReportWriter) {
}

ReportWriters::~ReportWriters() {
}

ScalarReportWriter* ReportWriters::getScalarReportWriter() const {
    return scalarReportWriter;
}


ArrayReportWriter* ReportWriters::getArrayBlockedReportWriter() const {
    return arrayReportWriter;
}

AccRejReportWriter* ReportWriters::getAccRejReportWriter() const {
    return accRejReportWriter;
}

