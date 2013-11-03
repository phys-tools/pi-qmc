#include "H5PartitionedScalarReportWriter.h"

#include <sstream>
#include "stats/ScalarEstimator.h"
#include "stats/PartitionedScalarAccumulator.h"
#include "stats/hdf5/H5Lib.h"

H5PartitionedScalarReportWriter::H5PartitionedScalarReportWriter(
        int nstep, hid_t writingGroupID)
:   nstep(nstep),
    writingGroupID(writingGroupID),
    groupID(0) {
}

H5PartitionedScalarReportWriter::~H5PartitionedScalarReportWriter() {
    delete [] groupID;
}

void H5PartitionedScalarReportWriter::createPartitionGroups() {
    groupID = new hid_t[partitionCount];
    for (int partition = 0; partition < partitionCount; ++partition) {
        std::string name = groupName(partition);
        groupID[partition] = H5Lib::createGroupInH5File(name, writingGroupID);
    }
}

void H5PartitionedScalarReportWriter::startReport(const ScalarEstimator *est,
        const PartitionedScalarAccumulator *acc) {
    partitionCount = acc->getPartitionCount();
    if (groupID == 0) createPartitionGroups();
    for (int partition = 0; partition < partitionCount; ++partition) {
        hid_t dataSetID
            = H5Lib::createScalarInH5File(*est, groupID[partition], nstep);
        datasetList.push_back(dataSetID);
    }
}

void H5PartitionedScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    datasetIterator = datasetList.begin();
}

void H5PartitionedScalarReportWriter::reportStep(const ScalarEstimator* est,
        const PartitionedScalarAccumulator *acc) {
    for (int partition = 0; partition < partitionCount; ++partition) {
        double value
            = (acc->getValue(partition) + est->getShift()) * est->getScale();
        H5Lib::writeScalarValue(*datasetIterator, istep, value);
        datasetIterator++;
    }
}

std::string H5PartitionedScalarReportWriter::groupName(int partition) {
    std::stringstream buffer;
    buffer << "partition" << partition;
    std::string name = std::string(buffer.str());
    return name;
}

