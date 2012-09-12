#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "H5ReportBuilder.h"
#include "stats/EstimatorIterator.h"
#include "stats/EstimatorManager.h"
#include "stats/ReportWriters.h"
#include "H5ArrayReportWriter.h"
#include "H5ScalarReportWriter.h"
#include "H5SplitScalarReportWriter.h"
#include "stats/NullAccRejReportWriter.h"

H5ReportBuilder::H5ReportBuilder(const std::string& filename,
        const EstimatorManager::SimInfoWriter *simInfoWriter)
:   filename(filename),
    simInfoWriter(simInfoWriter),
    fileID(0),
    writingGroupID(0),
    stepAttrID(0) {
    fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    simInfoWriter->writeH5(fileID);
}

H5ReportBuilder::~H5ReportBuilder() {
    delete reportWriters;
    H5Fclose(fileID);
}

void H5ReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->getNStep();
    istep = 0;
    writingGroupID = createH5Group("estimators", fileID);
    createReportWriters(manager);
    createStepAttribute();
    EstimatorIterator iterator = manager->getEstimatorIterator();
    do {
        (*iterator)->startReport(reportWriters);
    } while (iterator.step());
}

void H5ReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    scalarWriter->startBlock(istep);
    arrayWriter->startBlock(istep);
    EstimatorIterator iterator = manager->getEstimatorIterator();
    do {
        (*iterator)->reportStep(reportWriters);
    } while (iterator.step());
    H5Awrite(stepAttrID, H5T_NATIVE_INT, &(++istep));
    H5Fflush(fileID, H5F_SCOPE_LOCAL);
    if (istep == nstep) {
        closeDatasets();
    }
}

void H5ReportBuilder::recordInputDocument(const std::string &docstring) {
    hid_t dspaceID = H5Screate(H5S_SCALAR);
    hid_t dtypeID = H5Tcopy(H5T_C_S1);
    size_t size = docstring.size();
    H5Tset_size(dtypeID, size);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dsetID = H5Dcreate2(fileID, "simInfo/inputDoc", dtypeID, dspaceID,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    hid_t dsetID = H5Dcreate(fileID, "comment", dtypeID, dspaceID, H5P_DEFAULT);
#endif
    H5Dwrite(dsetID, dtypeID, H5S_ALL, H5S_ALL, H5P_DEFAULT, docstring.data());
    H5Dclose(dsetID);
    H5Tclose(dtypeID);
    H5Sclose(dspaceID);
}

hid_t H5ReportBuilder::createH5Group(std::string name, hid_t fileID) {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t groupID = H5Gcreate2(fileID, "estimators", H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
#else
    hid_t groupID = H5Gcreate(fileID,"estimators",0);
#endif 
    return groupID;
}

void H5ReportBuilder::createReportWriters(EstimatorManager*& manager) {
    scalarWriter = new H5ScalarReportWriter(nstep, writingGroupID);
    arrayWriter = new H5ArrayReportWriter(nstep, writingGroupID);
    NullAccRejReportWriter* accrejWriter = new NullAccRejReportWriter();
    reportWriters = new ReportWriters(scalarWriter, arrayWriter, accrejWriter);
}

void H5ReportBuilder::createStepAttribute() {
    hsize_t dims = 1;
    hid_t dataspaceID = H5Screate_simple(1, &dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    stepAttrID = H5Acreate2(writingGroupID, "nstep", H5T_NATIVE_INT,
            dataspaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
    stepAttrID = H5Acreate(writingGroupID, "nstep", H5T_NATIVE_INT,
            dataspaceID, H5P_DEFAULT);
#endif
    H5Sclose(dataspaceID);
    H5Awrite(stepAttrID, H5T_NATIVE_INT, &istep);
}

void H5ReportBuilder::closeDatasets() {
    H5Aclose (stepAttrID);
    H5Gclose (writingGroupID);
}

