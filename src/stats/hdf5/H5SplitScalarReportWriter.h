#ifndef H5SPLITSCALARREPORTWRITER_H_
#define H5SPLITSCALARREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
class ScalarEstimator;
#include <hdf5.h>
#include <vector>

class H5SplitScalarReportWriter
:   public ReportWriterInterface<ScalarEstimator> {
public:
    H5SplitScalarReportWriter(int nstep, hid_t writingGroupID);
    virtual ~H5SplitScalarReportWriter();

    virtual void startReport(const ScalarEstimator& est);
    virtual void startBlock(int istep);
    virtual void reportStep(const ScalarEstimator& est);

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
