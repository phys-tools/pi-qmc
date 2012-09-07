#ifndef __H5ReportBuilder_h_
#define __H5ReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include "EstimatorManager.h"
#include <string>
#include <vector>
#include <hdf5.h>

class EstimatorManager;
class ReportWriters;
class H5ScalarReportWriter;
class H5ArrayReportWriter;

/** Class for recording statistical data to an HDF5 file.

 The H5ReportBuilder is responsible for recording data
 recording it to disk. We have chosen an HDF5
 file format (http://hdf.ncsa.uiuc.edu/HDF5) so that we can compress
 the data and include meta data and other structure.
 @author John Shumway */
class H5ReportBuilder: public EstimatorReportBuilder {
public:
    H5ReportBuilder(const std::string &filename,
            const EstimatorManager::SimInfoWriter*);
    virtual ~H5ReportBuilder();

    virtual void initializeReport(EstimatorManager*);
    virtual void collectAndWriteDataBlock(EstimatorManager*);
    virtual void recordInputDocument(const std::string &docstring);

    virtual void startScalarReport(const ScalarEstimator& est);
    virtual void reportScalarStep(const ScalarEstimator& est);
    virtual void startArrayReport(const ArrayEstimator& est);
    virtual void reportArrayStep(const ArrayEstimator& est);

private:
    std::string filename;
    const EstimatorManager::SimInfoWriter *simInfoWriter;
    int nstep;
    int istep;
    typedef std::vector<hid_t> DataSetContainer;
    typedef DataSetContainer::iterator DataSetIter;

    hid_t fileID;
    hid_t writingGroupID;
    hid_t stepAttrID;

    DataSetContainer dataset;
    DataSetIter dset;

    ReportWriters *reportWriters;
    H5ScalarReportWriter *scalarWriter;
    H5ArrayReportWriter *arrayWriter;
};
#endif
