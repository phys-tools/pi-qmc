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
    initializeGrid(100);
    propagate();
    double value = readValue(100);
    return value;
}

void Propagator::setupGrid() {
    grid = new PropagatorGrid(256, 0.124235);
}

void Propagator::initializeGrid(int index0) {
    grid->initialize(index0);
}

void Propagator::propagate() {
    double deltaTau = 0.124235;
    double mass = 1.0;
    grid->setupKineticPropagator(mass, deltaTau);
    grid->setupPotentialPropagator(potential, deltaTau);
//    grid->evolveVHalfDeltaTau();
//    for (int i = 0; i < 100 - 2; ++i) {
//        grid->toKSpace();
//        grid->evolveTDeltaTau();
//        grid->toRealSpace();
//        grid->evolveVDeltaTau();
//    }
    grid->toKSpace();
    grid->evolveTDeltaTau();
    grid->toRealSpace();
//    grid->evolveVHalfDeltaTau();
}

double Propagator::readValue(int index) const {
    double value = grid->readValue(index);
    return value;
}

PropagatorGrid* Propagator::getGrid() const {
    return grid;
}

double Propagator::potential(double x) {
    return 0.0; //0.5 * x * x;
}




