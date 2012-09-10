#ifndef H5LIB_H_
#define H5LIB_H_
#include <hdf5.h>
#include <string>
class ScalarEstimator;

class H5Lib {
public:
    static hid_t createScalarDataset(int nstep, hid_t writingGroupID, std::string name);
    static void writeScalarValue(hid_t dataSetID, int istep, float value);
    static void writeUnitName(const std::string& unitName, hid_t& dataSetID);
    static void writeTypeString(const std::string& typeString, hid_t& dataSetID);
    static hid_t createScalarInH5File(const ScalarEstimator& est,
            hid_t writingGroupID, int nstep);
};

#endif
