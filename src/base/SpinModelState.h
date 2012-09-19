#ifndef __SpinModelState_h_
#define __SpinModelState_h_

#include "ModelState.h"
#include <cstdlib>
#include <blitz/array.h>

///Model state for spin up and spin down.
/// @version $Revision: 338 $
/// @author John Shumway and Jianheng Liu
class SpinModelState: public ModelState {
public:
    typedef blitz::Array<int, 1> IArray;
    SpinModelState(int npart, int initial);
    virtual ~SpinModelState() {
    }
    ;
    virtual void write(std::ostream &os) const;
    virtual bool read(const std::string &line);
    const IArray& getSpinState() const {
        return spinState;
    }
    IArray& getModelState() {
        return spinState;
    }
    virtual int getModelCount() const {
        return npart + 1;
    }
    virtual int getModelState() const;
    void flipSpin(int i) {
        spinState(i) = 1 - spinState(i);
    }
    virtual bool isSpinModelState() {
        return true;
    }
    virtual int getPartitionCount() const {
        return 1;
    }
    virtual double getValue(int i) const {
        return 1.0;
    }
    virtual void broadcastToMPIWorkers(const MPIManager *mpi);
private:
    const int npart;
    IArray spinState;
};
#endif
