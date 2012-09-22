#ifndef __EnumeratedModelState_h_
#define __EnumeratedModelState_h_

#include "ModelState.h"

class EnumeratedModelState: public ModelState {
public:
    EnumeratedModelState(int modelCount);
    virtual ~EnumeratedModelState() {
    }
    ;
    virtual void write(std::ostream &os) const;
    virtual bool read(const std::string &line);
    void setModelState(int i) {
        modelState = i;
    }
    virtual int getModelCount() const {
        return modelCount;
    }
    virtual int getModelState() const {
        return modelState;
    }
    virtual int getPartitionCount() const {
        return modelCount;
    }
    virtual void evaluate(Paths* paths) {
    }
    virtual double getValue(int i) const {
        return (i == modelState) ? 1.0 : 0.0;
    }
    virtual void broadcastToMPIWorkers(const MPIManager *mpi);
private:
    const int modelCount;
    int modelState;
};
#endif
