#ifndef __H5ReportBuilder_h_
#define __H5ReportBuilder_h_

#include "stats/EstimatorReportBuilder.h"
#include "stats/EstimatorManager.h"
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

private:
    std::string filename;
    const EstimatorManager::SimInfoWriter *simInfoWriter;
    int nstep;
    int istep;

    hid_t fileID;
    hid_t writingGroupID;
    hid_t stepAttrID;

    ReportWriters *reportWriters;
    H5ScalarReportWriter *scalarWriter;
    H5ArrayReportWriter *arrayWriter;
    void closeDatasets();
};
#endif
