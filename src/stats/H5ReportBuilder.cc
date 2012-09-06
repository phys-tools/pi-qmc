#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "H5ReportBuilder.h"
#include "ScalarEstimator.h"
#include "EstimatorManager.h"
#include "ArrayEstimator.h"
#include "ReportWriters.h"
#include "H5ScalarReportWriter.h"
#include "ReportWriterInterface.h"

H5ReportBuilder::H5ReportBuilder(const std::string& filename,
        const EstimatorManager::SimInfoWriter *simInfoWriter) :
        filename(filename), simInfoWriter(simInfoWriter), fileID(0), writingGroupID(
                0), stepAttrID(0), dataset(0) {
    fileID = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    simInfoWriter->writeH5(fileID);
    scalarWriter = new H5ScalarReportWriter();
//    arrayWriter = new H5ArrayReportWriter();
    ReportWriterInterface<AccRejEstimator> *accrejWriter
        = new ReportWriterInterface<AccRejEstimator>;
    reportWriters = new ReportWriters(scalarWriter, 0, 0);
}

H5ReportBuilder::~H5ReportBuilder() {
    H5Fclose(fileID);
}

void H5ReportBuilder::initializeReport(EstimatorManager *manager) {
    nstep = manager->nstep;
    scalarWriter->setNstep(nstep);
    istep = 0;
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    writingGroupID = H5Gcreate2(fileID, "estimators", H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
#else
    writingGroupID = H5Gcreate(fileID,"estimators",0);
#endif
    scalarWriter->setWritingGroupID(writingGroupID);
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
    scalarWriter->startScalarReport(est);
}

void H5ReportBuilder::reportScalarStep(const ScalarEstimator &est) {
    scalarWriter->reportScalarStep(est);
}

void H5ReportBuilder::startArrayReport(
        const ArrayEstimator& est) {
    hsize_t *dims;
    dims = new hsize_t[est.getNDim()];
    unsigned int maxDim = 1, imaxDim = 0, size = 1;
    // Find maximum dimension for compression.
    for (int i = 0; i < est.getNDim(); ++i) {
        dims[i] = est.getExtent(i);
        size *= dims[i];
        if (dims[i] > maxDim) {
            maxDim = dims[i];
            imaxDim = i;
        }
    }
    hid_t dataSpaceID = H5Screate_simple(est.getNDim(), dims, NULL);
    bool useCompression = (size > 10000);
    hid_t plist = H5P_DEFAULT;
    if (useCompression) {
        plist = H5Pcreate(H5P_DATASET_CREATE);
        dims[imaxDim] = dims[imaxDim] * 10000 / size;
        if (dims[imaxDim] == 0)
            dims[imaxDim] = 1;
        H5Pset_chunk(plist, est.getNDim(), dims);
        H5Pset_deflate(plist, 1);
    }
    delete dims;
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dataSetID = H5Dcreate2(writingGroupID, est.getName().c_str(),
            H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT, plist, H5P_DEFAULT);
#else
    hid_t dataSetID = H5Dcreate(writingGroupID, est.getName().c_str(),
            H5T_NATIVE_FLOAT, dataSpaceID, plist);
#endif
    {
        const std::string& typeString(est.getTypeString());
        hsize_t dims = 1;
        hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);
        hid_t strType = H5Tcopy(H5T_C_S1);
        H5Tset_size(strType, typeString.length());
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t attrID = H5Acreate2(dataSetID, "type", strType, dataSpaceID,
                H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t attrID = H5Acreate(dataSetID, "type", strType, dataSpaceID,
                H5P_DEFAULT);
#endif
        H5Awrite(attrID, strType, typeString.c_str());
        H5Sclose(dataSpaceID);
        H5Aclose(attrID);
    }
    dataset.push_back(dataSetID);
    if (est.hasError()) {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t dataSetID = H5Dcreate2(writingGroupID,
                (est.getName() + "_err").c_str(), H5T_NATIVE_FLOAT, dataSpaceID,
                H5P_DEFAULT, plist, H5P_DEFAULT);
#else
        hid_t dataSetID = H5Dcreate(writingGroupID, (est.getName()+"_err").c_str(),
                H5T_NATIVE_FLOAT, dataSpaceID, plist);
#endif
        dataset.push_back(dataSetID);
    }
    H5Sclose(dataSpaceID);
    if (useCompression)
        H5Pclose(plist);
    if (est.hasScale()) {
        hsize_t dim = est.getNDim();
        hid_t attrSpaceID = H5Screate_simple(1, &dim, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t attrID = H5Acreate2(dataSetID, "scale", H5T_NATIVE_DOUBLE,
                attrSpaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t attrID = H5Acreate(dataSetID, "scale",
                H5T_NATIVE_DOUBLE, attrSpaceID, H5P_DEFAULT);
#endif
        H5Sclose(attrSpaceID);
        H5Awrite(attrID, H5T_NATIVE_DOUBLE, est.getScale());
        H5Aclose(attrID);
    }

    if (est.hasOrigin()) {
        hsize_t dim = est.getNDim();
        hid_t attrSpaceID = H5Screate_simple(1, &dim, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t attrID = H5Acreate2(dataSetID, "origin", H5T_NATIVE_DOUBLE,
                attrSpaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t attrID = H5Acreate(dataSetID, "origin",
                H5T_NATIVE_DOUBLE, attrSpaceID, H5P_DEFAULT);
#endif
        H5Sclose(attrSpaceID);
        H5Awrite(attrID, H5T_NATIVE_DOUBLE, est.getOrigin());
        H5Aclose(attrID);
    }

    if (est.hasMin()) {
        hsize_t dim = est.getNDim();
        hid_t attrSpaceID = H5Screate_simple(1, &dim, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t attrID = H5Acreate2(dataSetID, "min", H5T_NATIVE_DOUBLE,
                attrSpaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t attrID = H5Acreate(dataSetID, "min",
                H5T_NATIVE_DOUBLE, attrSpaceID, H5P_DEFAULT);
#endif
        H5Sclose(attrSpaceID);
        H5Awrite(attrID, H5T_NATIVE_DOUBLE, est.getMin());
        H5Aclose(attrID);
    }

    if (est.hasMax()) {
        hsize_t dim = est.getNDim();
        hid_t attrSpaceID = H5Screate_simple(1, &dim, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t attrID = H5Acreate2(dataSetID, "max", H5T_NATIVE_DOUBLE,
                attrSpaceID, H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t attrID = H5Acreate(dataSetID, "max",
                H5T_NATIVE_DOUBLE, attrSpaceID, H5P_DEFAULT);
#endif
        H5Sclose(attrSpaceID);
        H5Awrite(attrID, H5T_NATIVE_DOUBLE, est.getMax());
        H5Aclose(attrID);
    }

}

void H5ReportBuilder::reportArrayStep(const ArrayEstimator& est) {
    est.normalize();
    H5Dwrite(*dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            est.getData());
    dset++;
    if (est.hasError()) {
        H5Dwrite(*dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                est.getError());
        dset++;
    }
    est.unnormalize();
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
