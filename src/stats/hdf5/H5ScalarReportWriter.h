#ifndef H5SCALARREPORTWRITER_H_
#define H5SCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
class ScalarEstimator;
class ScalarAccumulator;
#include <hdf5.h>
#include <vector>

class H5ScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator, ScalarAccumulator> {
public:
    H5ScalarReportWriter(int nstep, hid_t writingGroupID);
    virtual ~H5ScalarReportWriter();

    virtual void startReport(const ScalarEstimator* est,
            const ScalarAccumulator *acc);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator* est,
            const ScalarAccumulator *acc);

private:
    int nstep;
    int istep;
    hid_t writingGroupID;

    typedef std::vector<hid_t> DataSetContainer;
    typedef DataSetContainer::iterator DataSetIterator;
    DataSetContainer datasetList;
    DataSetIterator datasetIterator;
};

#endif
