#include "StdoutAccRejReportWriter.h"

StdoutAccRejReportWriter::StdoutAccRejReportWriter() {
}

StdoutAccRejReportWriter::~StdoutAccRejReportWriter() {
}

void StdoutAccRejReportWriter::reportStep(const AccRejEstimator& est) {
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
