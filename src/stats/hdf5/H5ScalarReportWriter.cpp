#include "H5ScalarReportWriter.h"

#include "stats/ScalarEstimator.h"
#include "stats/hdf5/H5Lib.h"

H5ScalarReportWriter::H5ScalarReportWriter()
:   nstep(0) {
}

H5ScalarReportWriter::~H5ScalarReportWriter() {
}

void H5ScalarReportWriter::startReport(const ScalarEstimator& est) {
    hid_t dataSetID = H5Lib::createScalarInH5File(est, writingGroupID, nstep);
    datasetList.push_back(dataSetID);
}

void H5ScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    datasetIterator = datasetList.begin();
}

void H5ScalarReportWriter::reportStep(const ScalarEstimator& est) {
    H5Lib::writeScalarValue(*datasetIterator, istep, est.getValue());
    datasetIterator++;
}

void H5ScalarReportWriter::setNstep(int nstep) {
    this->nstep = nstep;
}

void H5ScalarReportWriter::setWritingGroupID(hid_t writingGroupID) {
    this->writingGroupID = writingGroupID;
}
