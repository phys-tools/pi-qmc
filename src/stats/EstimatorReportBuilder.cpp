#include "EstimatorReportBuilder.h"

EstimatorReportBuilder::EstimatorReportBuilder() {
}

EstimatorReportBuilder::~EstimatorReportBuilder() {
}

void EstimatorReportBuilder::startScalarReport(const ScalarEstimator& est) {
}

void EstimatorReportBuilder::startArrayBlockedReport(
        const ArrayBlockedEstimator& est) {
}

void EstimatorReportBuilder::reportScalarStep(const ScalarEstimator& est) {
}

void EstimatorReportBuilder::startAccRejReport(const AccRejEstimator& est) {
}

void EstimatorReportBuilder::reportAccRejStep(const AccRejEstimator& est) {
}

void EstimatorReportBuilder::startArrayReport(const ArrayEstimator& est) {
}

void EstimatorReportBuilder::reportArrayStep(const ArrayEstimator& est) {
}

void EstimatorReportBuilder::reportArrayBlockedStep(
        const ArrayBlockedEstimator& est) {
}

void EstimatorReportBuilder::recordInputDocument(const std::string& docstring) {
}

