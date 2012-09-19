#ifndef H5SPLITSCALARREPORTWRITER_H_
#define H5SPLITSCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
class ScalarEstimator;
class PartitionedScalarAccumulator;
#include <hdf5.h>
#include <vector>
#include <string>

class H5PartitionedScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, PartitionedScalarAccumulator> {
public:
    H5PartitionedScalarReportWriter(int nstep, hid_t writingGroupID);
    virtual ~H5PartitionedScalarReportWriter();

    virtual void startReport(const ScalarEstimator *est,
            const PartitionedScalarAccumulator *acc);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator *est,
            const PartitionedScalarAccumulator *acc);

private:
    int nstep;
    int istep;
    int partitionCount;
    hid_t writingGroupID;

    typedef std::vector<hid_t> DataSetContainer;
    typedef DataSetContainer::iterator DataSetIterator;
    DataSetContainer datasetList;
    DataSetIterator datasetIterator;
    static std::string groupName(int partition);
};

#endif
