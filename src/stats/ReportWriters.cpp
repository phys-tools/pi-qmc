#include "ReportWriters.h"
#include "ReportWriterInterface.h"

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

void ReportWriters::startScalarReport(ScalarEstimator& est) {
    scalarReportWriter->startReport(est);
}

void ReportWriters::reportScalarStep(ScalarEstimator& est) {
    scalarReportWriter->reportStep(est);
}

void ReportWriters::startAccRejReport(AccRejEstimator& est) {
    accRejReportWriter->startReport(est);
}

void ReportWriters::reportAccRejStep(AccRejEstimator& est) {
    accRejReportWriter->reportStep(est);
}

void ReportWriters::startArrayReport(ArrayEstimator& est) {
    arrayReportWriter->startReport(est);
}

void ReportWriters::reportArrayStep(ArrayEstimator& est) {
    arrayReportWriter->reportStep(est);
}
