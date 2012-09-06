#ifndef H5SCALARREPORTWRITER_H_
#define H5SCALARREPORTWRITER_H_

#include "ReportWriterInterface.h"
class ScalarEstimator;
#include <hdf5.h>
#include <vector>

class H5ScalarReportWriter : public ReportWriterInterface<ScalarEstimator> {
public:
    H5ScalarReportWriter();
    virtual ~H5ScalarReportWriter();

    virtual void startScalarReport(const ScalarEstimator& est);
    virtual void reportScalarStep(const ScalarEstimator& est);

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
