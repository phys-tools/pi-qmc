#ifndef H5ARRAYREPORTWRITER_H_
#define H5ARRAYREPORTWRITER_H_

#include "stats/ReportWriterInterface.h"
#include "stats/ArrayEstimator.h"
#include <vector>
#include "hdf5.h"

class H5ArrayReportWriter
:   public ReportWriterInterface<ArrayEstimator, ScalarAccumulator> {
public:
    H5ArrayReportWriter(int nstep, hid_t writingGroupID);
    virtual ~H5ArrayReportWriter();

    virtual void startReport(const ArrayEstimator &est);
    virtual void startBlock(int istep);
    virtual void reportStep(const ArrayEstimator &est);

private:
    int nstep;
    int istep;
    hid_t writingGroupID;

    typedef std::vector<hid_t> DataSetContainer;
    typedef DataSetContainer::iterator DataSetIterator;
    DataSetContainer dataset;
    DataSetIterator dset;
};

#endif
