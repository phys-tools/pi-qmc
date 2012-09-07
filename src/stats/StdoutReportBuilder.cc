#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "StdoutReportBuilder.h"
#include "EstimatorManager.h"
#include "ReportWriters.h"
#include "StdoutAccRejReportWriter.h"
#include "StdoutArrayReportWriter.h"
#include "StdoutScalarReportWriter.h"
#include <ctime>
#include <string.h>
#include <iostream>

StdoutReportBuilder::StdoutReportBuilder()
:   istep(0),
    nstep(0) {
}

StdoutReportBuilder::~StdoutReportBuilder() {
}

void StdoutReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->nstep;
    int n = manager->estimator.size();

    scalarWriter = new StdoutScalarReportWriter(n, nstep);
    arrayWriter = new StdoutArrayReportWriter();
    accrejWriter = new StdoutAccRejReportWriter();
    reportWriters = new ReportWriters(scalarWriter, arrayWriter, accrejWriter);
}

void StdoutReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    writeBlockHeader();
    scalarWriter->startBlock(istep);
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->reportStep(*this);
    }
    std::cout << std::endl;
    ++istep;
}

void StdoutReportBuilder::writeBlockHeader() {
    std::time_t rawtime;
    std::time(&rawtime);
    char* myTime = std::ctime(&rawtime);
    myTime[strlen(myTime) - 1] = '\0'; //Hackish, gets rid of a trailing newline
    std::cout << "********** Block " << istep + 1 << " of " << nstep
            << " blocks **********";
    std::cout << " (" << myTime << ")" << std::endl;
}

void StdoutReportBuilder::reportScalarStep(const ScalarEstimator& est) {
    scalarWriter->reportStep(est);
}

void StdoutReportBuilder::reportAccRejStep(const AccRejEstimator& est) {
    accrejWriter->reportStep(est);
}

void StdoutReportBuilder::reportArrayBlockedStep(const ArrayEstimator& est) {
    arrayWriter->reportStep(est);
}
