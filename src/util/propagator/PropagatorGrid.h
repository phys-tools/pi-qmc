#ifndef PROPAGATORGRID_H_
#define PROPAGATORGRID_H_

class PropagatorGrid {
public:
    PropagatorGrid();
    virtual ~PropagatorGrid();

    void toRealSpace();
    void toKSpace();
    void evolveTDeltaTau();
    void evolveVDeltaTau();
    void evolveVHalfDeltaTau();
};

#endif
