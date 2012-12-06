#include "PropagatorGrid.h"

PropagatorGrid::PropagatorGrid(int size, double deltaX)
    :   size(size),
        deltaX(deltaX),
        deltaK(2.0 * PI / (deltaX * size)),
        value(new Complex[size]) {
}

PropagatorGrid::~PropagatorGrid() {
    delete value;
}

void PropagatorGrid::toRealSpace() {
}

void PropagatorGrid::toKSpace() {
}

void PropagatorGrid::evolveTDeltaTau() {
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




