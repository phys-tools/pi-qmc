#include "ReportWriters.h"

ReportWriters::ReportWriters(
        ReportWriterInterface<ScalarEstimator>* scalarReportWriter,
        ReportWriterInterface<ArrayEstimator>* arrayReportWriter,
        ReportWriterInterface<AccRejEstimator>* accRejReportWriter)
:   scalarReportWriter(scalarReportWriter),
    arrayReportWriter(arrayReportWriter),
    accRejReportWriter(accRejReportWriter) {
}

ReportWriters::~ReportWriters() {
    delete scalarReportWriter;
    delete arrayReportWriter;
    delete accRejReportWriter;
}

ReportWriterInterface<ScalarEstimator>* ReportWriters::getScalarReportWriter() const {
    return scalarReportWriter;
}


ReportWriterInterface<ArrayEstimator>* ReportWriters::getArrayBlockedReportWriter() const {
    return arrayReportWriter;
}

ReportWriterInterface<AccRejEstimator>* ReportWriters::getAccRejReportWriter() const {
    return accRejReportWriter;
}
