#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "StdoutReportBuilder.h"
#include "ScalarEstimator.h"
#include "AccRejEstimator.h"
#include "ArrayEstimator.h"
#include "EstimatorManager.h"
#include <ctime>
#include <string.h>
#include <iostream>

StdoutReportBuilder::StdoutReportBuilder() {
}

StdoutReportBuilder::~StdoutReportBuilder() {
}

void StdoutReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->nstep;
    istep = 0;
    int n = manager->estimator.size();
    sum.resize(n);
    sum2.resize(n);
    norm.resize(n);
    sum = 0;
    sum2 = 0;
    norm = 0;
}

void StdoutReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    std::time_t rawtime;
    std::time(&rawtime);
    char * myTime = std::ctime(&rawtime);
    myTime[strlen(myTime) - 1] = '\0'; //Hackish, gets rid of a trailing newline
    std::cout << "********** Block " << istep + 1 << " of " << nstep
            << " blocks **********"; //<< std::endl;
    std::cout << " (" << myTime << ")" << std::endl;
    iscalar = 0;
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->reportStep(*this);
    }
    std::cout << std::endl;
    if (++istep == nstep) {
    }
}

void StdoutReportBuilder::reportScalarStep(const ScalarEstimator& est) {
    double value = est.getValue();
    sum(iscalar) += value;
    sum2(iscalar) += value * value;
    norm(iscalar) += 1;
    std::cout << est.getName();
    if (est.getUnitName() != "")
        std::cout << " (" << est.getUnitName() << ")";
    std::cout << ": " << value << ", " << "Av="
            << sum(iscalar) / (norm(iscalar)) << " +-"
            << sqrt(
                    sum2(iscalar)
                            - sum(iscalar) * sum(iscalar) / (norm(iscalar)))
                    / (norm(iscalar) - 1) << std::endl;
    ++iscalar;
}

void StdoutReportBuilder::reportAccRejStep(const AccRejEstimator& est) {
    std::cout << est.getName() << std::endl;
    int nlevel = est.getNLevel();
    const AccRejEstimator::IArray& nacc(est.getNAccept());
    const AccRejEstimator::IArray& ntrl(est.getNTrial());
    for (int i = nlevel - 1; i >= 0; --i) {
        std::cout << "Level " << i << ": " << nacc(i) << "/" << ntrl(i) << " "
                << nacc(i) / (float) (ntrl(i) == 0 ? 1 : ntrl(i)) << std::endl;
    }
    std::cout << "Total:   " << nacc(0) << "/" << ntrl(nlevel - 1) << " "
            << nacc(0) / (float) (ntrl(nlevel - 1) == 0 ? 1 : ntrl(nlevel - 1))
            << std::endl;
}

void

StdoutReportBuilder::reportArrayBlockedStep(const ArrayEstimator& est) {
    std::cout << "(measured " << est.getName() << ")" << std::endl;
}
