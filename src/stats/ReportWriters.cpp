#include "ReportWriters.h"

ReportWriters::ReportWriters(ScalarReportWriter* scalarReportWriter,
        ArrayReportWriter* arrayReportWriter,
        ArrayBlockedReportWriter* arrayBlockedReportWriter,
        AccRejReportWriter* accRejReportWriter) :
        scalarReportWriter(scalarReportWriter), arrayReportWriter(
                arrayReportWriter), arrayBlockedReportWriter(
                arrayBlockedReportWriter), accRejReportWriter(
                accRejReportWriter) {
}

ReportWriters::~ReportWriters() {
}

ScalarReportWriter* ReportWriters::getScalarReportWriter() const {
    return scalarReportWriter;
}

ArrayReportWriter* ReportWriters::getArrayReportWriter() const {
    return arrayReportWriter;
}

ArrayBlockedReportWriter* ReportWriters::getArrayBlockedReportWriter() const {
    return arrayBlockedReportWriter;
}

AccRejReportWriter* ReportWriters::getAccRejReportWriter() const {
    return accRejReportWriter;
}

