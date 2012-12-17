#ifndef PROPAGATORGRID_H_
#define PROPAGATORGRID_H_

#include <complex>
class FFT1D;
class KineticGrid;
class PotentialGrid;

class PropagatorGrid {
public:
    typedef std::complex<double> Complex;

    PropagatorGrid(int size, double deltaX);
    virtual ~PropagatorGrid();

    void toRealSpace();
    void toKSpace();
    void evolveTDeltaTau();
    void evolveVDeltaTau();
    void evolveVHalfDeltaTau();

    double readValue() const;
    void initialize(int index0);
    Complex operator()(int index) const;

    void setupKineticPropagator(double mass, double deltaTau);
    void setupPotentialPropagator(double (*v)(double), double deltaTau);

    double getDeltaX() const;
    double getDeltaK() const;
private:
    const int size;
    double oneOverSqrtSize;
    const double deltaX;
    const double x0;
    const double deltaK;
    Complex* value;
    FFT1D* fft;
    KineticGrid* kineticPropagator;
    PotentialGrid* potentialPropagator;
    PotentialGrid* halfPotentialPropagator;

    static const double PI;
    void scaleBySqrtOfSize();
};

#endif
