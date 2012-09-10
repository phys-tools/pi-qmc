#include "H5SplitScalarReportWriter.h"

#include "stats/ScalarEstimator.h"
#include "stats/hdf5/H5Lib.h"

H5SplitScalarReportWriter::H5SplitScalarReportWriter(int nstep, hid_t writingGroupID)
:   nstep(nstep),
    writingGroupID(writingGroupID) {
}

H5SplitScalarReportWriter::~H5SplitScalarReportWriter() {
}

void H5SplitScalarReportWriter::startReport(const ScalarEstimator& est) {
    hid_t dataSetID = H5Lib::createScalarInH5File(est, writingGroupID, nstep);
    datasetList.push_back(dataSetID);
}

void H5SplitScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    datasetIterator = datasetList.begin();
}

void H5SplitScalarReportWriter::reportStep(const ScalarEstimator& est) {
    H5Lib::writeScalarValue(*datasetIterator, istep, est.getValue());
    datasetIterator++;
}
