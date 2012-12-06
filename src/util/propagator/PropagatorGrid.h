#ifndef PROPAGATORGRID_H_
#define PROPAGATORGRID_H_

#include <complex>

class PropagatorGrid {
public:
    typedef std::complex<double> Complex;

    PropagatorGrid(int size);
    virtual ~PropagatorGrid();

    void toRealSpace();
    void toKSpace();
    void evolveTDeltaTau();
    void evolveVDeltaTau();
    void evolveVHalfDeltaTau();

    double readValue() const;
    void initialize(int index0);
    Complex operator()(int index) const;
private:
    const int size;
    Complex* value;
};

#endif
