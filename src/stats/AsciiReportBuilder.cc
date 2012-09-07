#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "AsciiReportBuilder.h"
#include "EstimatorManager.h"
#include "ReportWriters.h"
#include "AsciiScalarReportWriter.h"
#include "NullArrayReportWriter.h"
#include "NullAccRejReportWriter.h"
#include <fstream>
#include <iostream>

AsciiReportBuilder::AsciiReportBuilder(const std::string& filename)
:   file(filename.c_str()) {
    scalarWriter = new AsciiScalarReportWriter(file);
    NullArrayReportWriter *arrayWriter = new NullArrayReportWriter();
    NullAccRejReportWriter *accrejWriter = new NullAccRejReportWriter();
    reportWriters = new ReportWriters(scalarWriter, arrayWriter, accrejWriter);
}

AsciiReportBuilder::~AsciiReportBuilder() {
    delete reportWriters;
}

void AsciiReportBuilder::initializeReport(EstimatorManager *manager) {
    file << "#";
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->startReport(reportWriters);
    }
    file << std::endl;
}

void AsciiReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->reportStep(reportWriters);
    }
    file << std::endl;
}
