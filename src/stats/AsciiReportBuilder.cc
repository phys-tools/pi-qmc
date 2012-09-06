#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "AsciiReportBuilder.h"
#include "ScalarEstimator.h"
#include "EstimatorManager.h"
#include <fstream>
#include <iostream>

AsciiReportBuilder::AsciiReportBuilder(const std::string& filename) :
        file(filename.c_str()) {
}

AsciiReportBuilder::~AsciiReportBuilder() {
}

void AsciiReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->nstep;
    istep = 0;
    file << "#";
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        const std::string typeString = (*est)->getTypeString();
        if (typeString.substr(0, 6) == "scalar") {
            file << "\"" << (*est)->getName();

            const std::string& unitName((*est)->getUnitName());
            if (unitName != "")
                file << " (" << unitName << ")";
            file << "\" ";
        }
    }
    file << std::endl;
}

void AsciiReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->reportStep(*this);
    }
    file << std::endl;
}

void AsciiReportBuilder::reportScalarStep(const ScalarEstimator& est) {
    file << est.getValue() << " ";
}
