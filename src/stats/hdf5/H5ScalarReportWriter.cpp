#include "H5ScalarReportWriter.h"

#include "stats/ScalarEstimator.h"
#include "stats/hdf5/H5Lib.h"

H5ScalarReportWriter::H5ScalarReportWriter(int nstep, hid_t writingGroupID)
:   nstep(nstep),
    writingGroupID(writingGroupID) {
}

H5ScalarReportWriter::~H5ScalarReportWriter() {
}

void H5ScalarReportWriter::startReport(const ScalarEstimator *est,
        const SimpleScalarAccumulator *acc) {
    hid_t dataSetID = H5Lib::createScalarInH5File(*est, writingGroupID, nstep);
    datasetList.push_back(dataSetID);
}

void H5ScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    datasetIterator = datasetList.begin();
}

void H5ScalarReportWriter::reportStep(const ScalarEstimator *est,
        const SimpleScalarAccumulator *acc) {
    H5Lib::writeScalarValue(*datasetIterator, istep, est->getValue());
    datasetIterator++;
}
