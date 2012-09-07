#ifndef H5ARRAYREPORTWRITER_H_
#define H5ARRAYREPORTWRITER_H_

#include "ReportWriterInterface.h"
#include "ArrayEstimator.h"
#include <vector>
#include "hdf5.h"

class H5ArrayReportWriter: public ReportWriterInterface<ArrayEstimator> {
public:
    H5ArrayReportWriter();
    virtual ~H5ArrayReportWriter();

    virtual void startReport(const ArrayEstimator &est);
    virtual void reportStep(const ArrayEstimator &est);

    void setNstep(int nstep);
    void setWritingGroupID(hid_t writingGroupID);
    void startBlock(int istep);

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
