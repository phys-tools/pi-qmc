#ifndef __SimInfoWriter_h_
#define __SimInfoWriter_h_

#include <string>
#include <iostream>
#include <vector>
#include "stats/EstimatorManager.h"
class SimulationInfo;

/// Class for writing out the SimulationInfo.
/// used in parsing and setup. 
class SimInfoWriter: public EstimatorManager::SimInfoWriter {
public:
    SimInfoWriter(const SimulationInfo &simInfo);
    virtual ~SimInfoWriter() {}
    virtual void writeH5(hid_t) const;
    //virtual void writeStdout(std::ostream&) const {};
    //virtual void writeAscii(std::ostream&) const {};
private:
    const SimulationInfo &simInfo;
};
#endif
