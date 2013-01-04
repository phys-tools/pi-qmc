#include "GridParameters.h"
#include <cmath>


GridParameters::GridParameters() {
    gridCount = 2048;
    deltaX = 0.005;
    xmin = -5.11;
    index0 = 1222;
}


GridParameters::GridParameters(double mass, double tau, double x0,
        double deltaX) {
    double width = 8.0 * calculateThermalWidth(mass, tau) / deltaX;
    gridCount = powerTwoCeiling(width);
    index0 = gridCount / 2;
    this->deltaX = deltaX;
    xmin = x0 - index0 * deltaX;
 }

GridParameters::~GridParameters() {
}

int GridParameters::getGridCount() const {
    return gridCount;
}

double GridParameters::getDeltaX() const {
    return deltaX;
}

double GridParameters::getXMin() const {
    return xmin;
}

int GridParameters::getIndex0() const {
    return index0;
}

double GridParameters::calculateThermalWidth(double mass, double tau) {
    return sqrt(tau / mass);
}
int GridParameters::powerTwoCeiling(double min) {
    int temp = 1;
    while (min > 1.0) {
        min /= 2.0;
        temp <<= 1;
    }
    return temp;
}

