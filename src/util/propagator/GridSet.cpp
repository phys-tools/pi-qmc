#include "GridSet.h"
#include "PropagatorGrid.h"
#include "GridParameters.h"

GridSet::GridSet() {

}

GridSet::~GridSet() {
    for (GridIterator iter = grid.begin(); iter != grid.end(); ++ iter) {
        delete *iter;
    }
    delete parameters;
}

double GridSet::getDeltaX() const {
    return parameters->getDeltaX();
}

void GridSet::setupGrid() {
    int gridCount = parameters->getGridCount();
    double deltaX = parameters->getDeltaX();
    double xmin = parameters->getXMin();
    PropagatorGrid* propagatorGrid = new PropagatorGrid(gridCount, deltaX,
            xmin);
    propagatorGrid->initialize(parameters->getIndex0());
    grid.push_back(propagatorGrid);
}

void GridSet::initializeGrid() {
    int index0 = parameters->getIndex0();
    for (GridIterator iter = grid.begin(); iter != grid.end(); ++ iter) {
        (*iter)->initialize(index0);
    }
}

void GridSet::setGridParameters(GridParameters* parameters) {
    this->parameters = parameters;
}

double GridSet::readValue0() const {
    return grid[grid.size() - 1]->readValue(parameters->getIndex0());
}

bool GridSet::isConverged() const {
    return grid.size() > 40;
}

PropagatorGrid* GridSet::getGrid(int index) {
    if (index >= grid.size()) {
        setupGrid();
        return getGrid(index);
    } else {
        return grid[index];
    }
}




