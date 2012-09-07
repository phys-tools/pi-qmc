#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "H5ReportBuilder.h"
#include "EstimatorManager.h"
#include "ReportWriters.h"
#include "H5ArrayReportWriter.h"
#include "H5ScalarReportWriter.h"
#include "NullAccRejReportWriter.h"

H5ReportBuilder::H5ReportBuilder(const std::string& filename,
        const EstimatorManager::SimInfoWriter *simInfoWriter) :
        filename(filename), simInfoWriter(simInfoWriter), fileID(0), writingGroupID(
                0), stepAttrID(0), dataset(0) {
    fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    simInfoWriter->writeH5(fileID);
    scalarWriter = new H5ScalarReportWriter();
    arrayWriter = new H5ArrayReportWriter();
    NullAccRejReportWriter *accrejWriter = new NullAccRejReportWriter();
    reportWriters = new ReportWriters(scalarWriter, 0, 0);
}

H5ReportBuilder::~H5ReportBuilder() {
    delete reportWriters;
    H5Fclose(fileID);
}

void H5ReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->nstep;
    scalarWriter->setNstep(nstep);
    arrayWriter->setNstep(nstep);
    istep = 0;
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    writingGroupID = H5Gcreate2(fileID, "estimators", H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
#else
    writingGroupID = H5Gcreate(fileID,"estimators",0);
#endif
    scalarWriter->setWritingGroupID(writingGroupID);
    arrayWriter->setWritingGroupID(writingGroupID);
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
    dataset.resize(0);
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->startReport(*this);
    }
}

void H5ReportBuilder::collectAndWriteDataBlock(EstimatorManager *manager) {
    dset = dataset.begin();
    scalarWriter->startBlock(istep);
    arrayWriter->startBlock(istep);
    for (EstimatorManager::EstimatorIter est = manager->estimator.begin();
            est != manager->estimator.end(); ++est) {
        (*est)->reportStep(*this);
    }
    H5Awrite(stepAttrID, H5T_NATIVE_INT, &(++istep));
    H5Fflush(fileID, H5F_SCOPE_LOCAL);
    if (istep == nstep) {
        H5Aclose(stepAttrID);
        for (DataSetIter d = dataset.begin(); d != dataset.end(); ++d) {
            H5Dclose(*d);
        }
        H5Gclose(writingGroupID);
        dataset.resize(0);
    }
}

void H5ReportBuilder::startScalarReport(const ScalarEstimator &est) {
    scalarWriter->startReport(est);
}

void H5ReportBuilder::reportScalarStep(const ScalarEstimator &est) {
    scalarWriter->reportStep(est);
}

void H5ReportBuilder::startArrayReport(const ArrayEstimator& est) {
    arrayWriter->startReport(est);
}

void H5ReportBuilder::reportArrayStep(const ArrayEstimator& est) {
    arrayWriter->reportStep(est);
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
