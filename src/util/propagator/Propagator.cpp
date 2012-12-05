#include "Propagator.h"
#include "PropagatorGrid.h"

Propagator::Propagator()
    :   grid(0) {
}

Propagator::~Propagator() {
    delete grid;
}

double Propagator::evaluate() {
    setupGrid();
    initializeGrid();
    propagate();
    double value = readValue();
    return value;
}

void Propagator::setupGrid() {
    grid = new PropagatorGrid();
}

void Propagator::initializeGrid() {
}

void Propagator::propagate() {
}

double Propagator::readValue() const {
    double value = grid->readValue();
    return value;
}

PropagatorGrid* Propagator::getGrid() const {
    return grid;
}



