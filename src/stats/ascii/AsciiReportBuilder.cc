#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "AsciiReportBuilder.h"
#include "stats/EstimatorManager.h"
#include "stats/ReportWriters.h"
#include "stats/ascii/AsciiScalarReportWriter.h"
#include "stats/NullArrayReportWriter.h"
#include "stats/NullAccRejReportWriter.h"
#include "stats/EstimatorIterator.h"
#include <fstream>

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
    EstimatorIterator iterator = manager->getEstimatorIterator();
    do {
        (*iterator)->startReport(reportWriters);
    } while (iterator.step());
    file << std::endl;
}

void AsciiReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    EstimatorIterator iterator = manager->getEstimatorIterator();
    do {
        (*iterator)->reportStep(reportWriters);
    } while (iterator.step());
    file << std::endl;
}
