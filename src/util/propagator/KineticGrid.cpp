#include "KineticGrid.h"
#include <cmath>

KineticGrid::KineticGrid(int gridCount, double deltaK, double mass,
        double deltaTau)
    :   gridCount(gridCount),
        value(new double[gridCount]) {
    value[0] = 1.0;
    for (int index = 1; index <= gridCount / 2; ++index) {
        double k = deltaK * index;
        double energy = 0.5 * k * k / mass;
        value[index] = value[gridCount - index] = exp(-deltaTau * energy);
    }
}

KineticGrid::~KineticGrid() {
    delete [] value;
}

double KineticGrid::operator ()(int index) const {
    return value[index];
}


