#include "PropagatorGrid.h"
#include "util/fft/FFT1D.h"
#include "util/propagator/KineticGrid.h"
#include "util/propagator/PotentialGrid.h"

PropagatorGrid::PropagatorGrid(int size, double deltaX)
    :   size(size),
        oneOverSqrtSize(1.0 / sqrt(size)),
        deltaX(deltaX),
        x0(0.0),
        deltaK(2.0 * PI / (deltaX * size)),
        value(new Complex[size]),
        fft(new FFT1D(value, size)),
        kineticPropagator(0),
        potentialPropagator(0),
        halfPotentialPropagator(0) {
}

PropagatorGrid::~PropagatorGrid() {
    delete fft;
    delete value;
    delete kineticPropagator;
    delete potentialPropagator;
    delete halfPotentialPropagator;
}

void PropagatorGrid::toRealSpace() {
    fft->reverse();
    scaleBySqrtOfSize();
}

void PropagatorGrid::toKSpace() {
    fft->forward();
    scaleBySqrtOfSize();
}

void PropagatorGrid::setupKineticPropagator(double mass, double deltaTau) {
    kineticPropagator = new KineticGrid(size, deltaK, mass, deltaTau);
}

void PropagatorGrid::setupPotentialPropagator(double (*v)(double),
        double deltaTau) {
    potentialPropagator =
            new PotentialGrid(size, deltaX, x0, v, deltaTau);
    halfPotentialPropagator =
            new PotentialGrid(size, deltaX, x0, v, 0.5 * deltaTau);
}

void PropagatorGrid::scaleBySqrtOfSize() {
    for (int i = 0; i < size; ++i) {
        value[i] *= oneOverSqrtSize;
    }
}

void PropagatorGrid::evolveTDeltaTau() {
    for (int i = 0; i < size; ++i) {
        value[i] *= (*kineticPropagator)(i);
    }
}

void PropagatorGrid::evolveVDeltaTau() {
    for (int i = 0; i < size; ++i) {
        value[i] *= (*potentialPropagator)(i);
    }
}

void PropagatorGrid::evolveVHalfDeltaTau() {
    for (int i = 0; i < size; ++i) {
        value[i] *= (*halfPotentialPropagator)(i);
    }
}

double PropagatorGrid::readValue(int index) const {
//    return 0.23182433490559262;
    return real(value[index]);
}

void PropagatorGrid::initialize(int index0) {
    for (int i = 0; i < size; ++i) {
        value[i] = 0.0;
    }
    value[index0] = 1.0;
}

PropagatorGrid::Complex PropagatorGrid::operator()(int index) const {
    return value[index];
}

const double PropagatorGrid::PI = 3.141592653589793;

double PropagatorGrid::getDeltaX() const {
    return deltaX;
}

double PropagatorGrid::getDeltaK() const {
    return deltaK;
}

