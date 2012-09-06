#ifndef __H5ReportBuilder_h_
#define __H5ReportBuilder_h_

#include "EstimatorReportBuilder.h"
#include <string>
#include <vector>
#include <hdf5.h>

#include "EstimatorManager.h"
class ReportWriters;
class H5ScalarReportWriter;

/** Class for recording statistical data to an HDF5 file.

 The H5ReportBuilder s responsible for recording data
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

    /// Method to start reporting a ScalarEstimator.
    virtual void startScalarReport(const ScalarEstimator& est);
    /// Method to write a step for a ScalarEstimator.
    virtual void reportScalarStep(const ScalarEstimator& est);
    /// Method to start reporting a ArrayBlockedEstimator.
    virtual void startArrayBlockedReport(const ArrayBlockedEstimator& est);
    /// Method to write a step for a ArrayBlockedEstimator.
    virtual void reportArrayBlockedStep(const ArrayBlockedEstimator& est);
    /// Method to start reporting a AccRejEstimator.
    virtual void startAccRejReport(const AccRejEstimator& est);
    /// Method to write a step for a AccRejEstimator.
    virtual void reportAccRejStep(const AccRejEstimator& est);
    /// Method to record the input file contents.
    virtual void recordInputDocument(const std::string &docstring);
private:
    std::string filename;
    const EstimatorManager::SimInfoWriter *simInfoWriter;
    int nstep;
    int istep;
    typedef std::vector<hid_t> DataSetContainer;
    typedef DataSetContainer::iterator DataSetIter;
    /// The HDF5 file.
    hid_t fileID;
    /// The current HDF5 writing group.
    hid_t writingGroupID;
    /// The group's step counter.
    hid_t stepAttrID;
    /// The current HDF5 datasets.
    DataSetContainer dataset;
    /// Iterator pointing at the current dataset.
    DataSetIter dset;
    ReportWriters *reportWriters;
    H5ScalarReportWriter *scalarWriter;
};
#endif
