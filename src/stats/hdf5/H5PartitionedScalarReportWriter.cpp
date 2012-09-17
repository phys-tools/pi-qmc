#include "H5PartitionedScalarReportWriter.h"

#include "stats/ScalarEstimator.h"
#include "stats/PartitionedScalarAccumulator.h"
#include "stats/hdf5/H5Lib.h"

H5PartitionedScalarReportWriter::H5PartitionedScalarReportWriter(
        int nstep, hid_t writingGroupID)
:   nstep(nstep),
    writingGroupID(writingGroupID) {
}

H5PartitionedScalarReportWriter::~H5PartitionedScalarReportWriter() {
}

void H5PartitionedScalarReportWriter::startReport(const ScalarEstimator *est,
        const PartitionedScalarAccumulator *acc) {
    hid_t dataSetID = H5Lib::createScalarInH5File(*est, writingGroupID, nstep);
    datasetList.push_back(dataSetID);
}

void H5PartitionedScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    datasetIterator = datasetList.begin();
}

void H5PartitionedScalarReportWriter::reportStep(const ScalarEstimator* est,
        const PartitionedScalarAccumulator *acc) {
    H5Lib::writeScalarValue(*datasetIterator, istep, est->getValue());
    datasetIterator++;
}
