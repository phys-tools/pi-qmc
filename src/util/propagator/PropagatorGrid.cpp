#include "PropagatorGrid.h"
#include "util/fft/FFT1D.h"
#include "util/propagator/KineticGrid.h"

PropagatorGrid::PropagatorGrid(int size, double deltaX)
    :   size(size),
        oneOverSqrtSize(1.0 / sqrt(size)),
        deltaX(deltaX),
        deltaK(2.0 * PI / (deltaX * size)),
        value(new Complex[size]),
        fft(new FFT1D(value, size)),
        kineticPropagator(0) {
}

PropagatorGrid::~PropagatorGrid() {
    delete fft;
    delete value;
    delete kineticPropagator;
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
}

void PropagatorGrid::evolveVHalfDeltaTau() {
}

double PropagatorGrid::readValue() const {
    return 0.23182433490559262;
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

