#include "H5ArrayReportWriter.h"

H5ArrayReportWriter::H5ArrayReportWriter()
:   nstep(nstep) {
}

H5ArrayReportWriter::~H5ArrayReportWriter() {
}

void H5ArrayReportWriter::startReport(const ArrayEstimator& est) {
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

void H5ArrayReportWriter::reportStep(const ArrayEstimator& est) {
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

void H5ArrayReportWriter::setWritingGroupID(hid_t writingGroupID) {
    this->writingGroupID = writingGroupID;
}

void H5ArrayReportWriter::setNstep(int nstep) {
    this->nstep = nstep;
}

void H5ArrayReportWriter::startBlock(int istep) {
    this->istep = istep;
    dset = dataset.begin();
}
