#include "PotentialGrid.h"
#include <cmath>

PotentialGrid::PotentialGrid(int size, double deltaX, double x0,
        double (*v)(double), double deltaTau)
    :   value(new double[size]),
        size(size) {
    for (int i = 0; i < size; ++i) {
        double x = x0 + i * deltaX;
        value[i] = exp(-deltaTau * v(x));
    }
}

PotentialGrid::~PotentialGrid() {
    delete value;
}

double PotentialGrid::operator()(int index) const {
    return value[index];
}

