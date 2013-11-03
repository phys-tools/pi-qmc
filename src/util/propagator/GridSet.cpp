#include "GridSet.h"
#include "PropagatorGrid.h"
#include "GridParameters.h"
#include "util/math/VPolyFit.h"

GridSet::GridSet()
    :   gridCount(0),
        lastDelta(0.0) {
}

GridSet::~GridSet() {
    for (int index = 0; index < gridCount; ++index) {
        delete gridArray[index];
    }
    delete parameters;
}

double GridSet::getDeltaX() const {
    return parameters->getDeltaX();
}

PropagatorGrid* GridSet::allocateGrid() {
    int gridCount = parameters->getGridCount();
    double deltaX = parameters->getDeltaX();
    double xmin = parameters->getXMin();
    PropagatorGrid* propagatorGrid = new PropagatorGrid(gridCount, deltaX,
            xmin);
    propagatorGrid->initialize(parameters->getIndex0());
    return propagatorGrid;
}


void GridSet::setGridParameters(GridParameters* parameters) {
    this->parameters = parameters;
}

void GridSet::extrapolateValue0() {
    int index0 = parameters->getIndex0();
    Array value(gridArray.size());
    for (int index = 0; index < gridCount; ++index) {
        value[index] = gridArray[index]->readValue(index0);
    }
    VPolyFit fitter(gridCount, 1, &deltaTau[0], &value[0]);
    fitter.fit();
    lastDelta = *fitter.getLastDelta();
    solution = *fitter.getSolution();
}

double GridSet::readValue0() const{
    return solution;
}

bool GridSet::isConverged() const {
    return fabs(lastDelta) < TOLERANCE;
}

PropagatorGrid* GridSet::makeNewGrid(double deltaTau) {
    PropagatorGrid* newGrid = allocateGrid();
    gridArray.push_back(newGrid);
    ++gridCount;
    this->deltaTau.push_back(deltaTau);
    return newGrid;
}

const double GridSet::TOLERANCE = 1e-12;
