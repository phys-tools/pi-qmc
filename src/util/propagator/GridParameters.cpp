#include "GridParameters.h"

GridParameters::GridParameters() {
    gridCount = 2048;
    deltaX = 0.005;
    xmin = -5.11;
    index0 = 1224;
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
