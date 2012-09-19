#include "H5Lib.h"
#include "stats/ScalarEstimator.h"

hid_t H5Lib::createScalarInH5File(const ScalarEstimator& est,
        hid_t writingGroupID, int nstep) {
    hid_t dataSetID = createScalarDataset(nstep, writingGroupID, est.getName());
    writeTypeString(est.getTypeString(), dataSetID);
    writeUnitName(est.getUnitName(), dataSetID);
    return dataSetID;
}

hid_t H5Lib::createScalarDataset(int nstep, hid_t writingGroupID,
        std::string name) {
    hsize_t dims = nstep;
    hid_t dataSpaceID = H5Screate_simple(1, &dims, NULL);
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t dataSetID = H5Dcreate2(writingGroupID, name.c_str(),
        H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT);
#else
    hid_t dataSetID = H5Dcreate(writingGroupID, name.c_str(),
        H5T_NATIVE_FLOAT, dataSpaceID, H5P_DEFAULT);
#endif
    H5Sclose(dataSpaceID);
    return dataSetID;
}

void H5Lib::writeTypeString(const std::string& typeString,
        hid_t& dataSetID) {
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

void H5Lib::writeUnitName(const std::string& unitName,
        hid_t& dataSetID) {
    if (unitName == "") return;
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

void H5Lib::writeScalarValue(hid_t dataSetID, int istep,
        float value) {
    hid_t space = H5Dget_space(dataSetID);
    H5Sselect_none(space);
    hsize_t element = istep;
    H5Sselect_elements(space, H5S_SELECT_SET, 1, &element);
    hsize_t dims = 1, maxdims = 1;
    hid_t mspace = H5Screate_simple(1, &dims, &maxdims);
    H5Dwrite(dataSetID, H5T_NATIVE_FLOAT, mspace, space, H5P_DEFAULT, &value);
}

hid_t H5Lib::createGroupInH5File(const std::string& name,
        hid_t containerID) {
#if (H5_VERS_MAJOR>1)||((H5_VERS_MAJOR==1)&&(H5_VERS_MINOR>=8))
    hid_t groupID = H5Gcreate2(containerID, name.c_str(), H5P_DEFAULT, H5P_DEFAULT,
            H5P_DEFAULT);
#else
    hid_t groupID = H5Gcreate(containerID, name.c_str(), 0);
#endif
    return groupID;
}

