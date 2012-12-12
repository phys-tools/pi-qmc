#include "PotentialGrid.h"

PotentialGrid::PotentialGrid(int size, double deltaX, double x0,
        double (*v)(double))
    :   value(new double[size]),
        size(size) {
}

PotentialGrid::~PotentialGrid() {
    delete value;
}
