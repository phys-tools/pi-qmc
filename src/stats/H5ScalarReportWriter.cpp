#include "H5ScalarReportWriter.h"

#include "ScalarEstimator.h"

H5ScalarReportWriter::H5ScalarReportWriter() :
        nstep(0) {
}

H5ScalarReportWriter::~H5ScalarReportWriter() {
}

void H5ScalarReportWriter::startScalarReport(const ScalarEstimator& est) {
    hsize_t dims = nstep;
    hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dataSetID = H5Dcreate2(writingGroupID, est.getName().c_str(),
            H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
#else
    hid_t dataSetID = H5Dcreate(writingGroupID, est.getName().c_str(),
            H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT);
#endif
    H5Sclose(dataSpaceID);
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
    const std::string& unitName(est.getUnitName());
    if (unitName != "") {
        hsize_t dims = 1;
        hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);
        hid_t strType = H5Tcopy(H5T_C_S1);
        H5Tset_size(strType, unitName.length());
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
        hid_t attrID = H5Acreate2(dataSetID, "unit", strType, dataSpaceID,
                H5P_DEFAULT, H5P_DEFAULT);
#else
        hid_t attrID = H5Acreate(dataSetID, "unit", strType, dataSpaceID,
                H5P_DEFAULT);
#endif
        H5Awrite(attrID, strType, unitName.c_str());
        H5Sclose(dataSpaceID);
        H5Aclose(attrID);
    }
    dataset.push_back(dataSetID);
}

void H5ScalarReportWriter::reportScalarStep(const ScalarEstimator& est) {
    hid_t space = H5Dget_space(*dset);
    H5Sselect_none(space);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR==6)&&(H5_VERS_RELEASE>=7))
    hsize_t element = istep;
    H5Sselect_elements(space, H5S_SELECT_SET, 1, &element);
#else
    hsize_t element[][1] = { {istep}};
    H5Sselect_elements(space, H5S_SELECT_SET, 1, (const hsize_t **)element);
#endif

    float a = est.getValue();
    hsize_t dims = 1, maxdims = 1;
    hid_t mspace = H5Screate_simple(1, &dims, &maxdims);
    H5Dwrite(*dset, H5T_NATIVE_FLOAT, mspace, space, H5P_DEFAULT, &a);
    dset++;
}

void H5ScalarReportWriter::setWritingGroupID(hid_t writingGroupID) {
    this->writingGroupID = writingGroupID;
}

void H5ScalarReportWriter::setNstep(int nstep) {
    this->nstep = nstep;
}

void H5ScalarReportWriter::startBlock(int istep) {
    this->istep = istep;
    dset = dataset.begin();
}

