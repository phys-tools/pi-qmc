#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "StdoutReportBuilder.h"
#include "stats/EstimatorManager.h"
#include "stats/EstimatorIterator.h"
#include "stats/ReportWriters.h"
#include "StdoutAccRejReportWriter.h"
#include "StdoutArrayReportWriter.h"
#include "StdoutScalarReportWriter.h"
#include <ctime>
#include <string.h>
#include <iostream>

StdoutReportBuilder::StdoutReportBuilder()
:   istep(0),
    nstep(0),
    reportWriters(0) {
}

StdoutReportBuilder::~StdoutReportBuilder() {
    delete reportWriters;
}

void StdoutReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->getNStep();
    scalarWriter = new StdoutScalarReportWriter(nstep);
    arrayWriter = new StdoutArrayReportWriter();
    accrejWriter = new StdoutAccRejReportWriter();
    reportWriters = new ReportWriters(scalarWriter, arrayWriter, accrejWriter);
    EstimatorIterator iterator = manager->getEstimatorIterator();
    do {
        (*iterator)->startReport(reportWriters);
    } while(iterator.step());}

void StdoutReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    writeBlockHeader();
    scalarWriter->startBlock(istep);
    EstimatorIterator iterator = manager->getEstimatorIterator();
    do {
        (*iterator)->reportStep(reportWriters);
    } while(iterator.step());
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

