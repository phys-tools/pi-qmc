#include "Propagator.h"
#include "PropagatorGrid.h"

Propagator::Propagator()
    :   grid(0) {
    potential = zeroPotential;
}

Propagator::~Propagator() {
    delete grid;
}

double Propagator::evaluate() {
    setupGrid();
    initializeGrid(1224);
    propagate();
    double value = readValue(1224);
    return value;
}

void Propagator::setupGrid() {
    grid = new PropagatorGrid(2048, 0.005, -5.11);
}

void Propagator::initializeGrid(int index0) {
    grid->initialize(index0);
}

void Propagator::propagate() {
    int stepCount = 100;
    double deltaTau = 0.124235 / stepCount;
    double mass = 1.0;
    grid->setupKineticPropagator(mass, deltaTau);
    grid->setupPotentialPropagator(potential, deltaTau);
    grid->evolveVHalfDeltaTau();
    for (int i = 0; i < stepCount - 1; ++i) {
        grid->toKSpace();
        grid->evolveTDeltaTau();
        grid->toRealSpace();
        grid->evolveVDeltaTau();
    }
    grid->toKSpace();
    grid->evolveTDeltaTau();
    grid->toRealSpace();
    grid->evolveVHalfDeltaTau();
}

double Propagator::readValue(int index) const {
    double value = grid->readValue(index);
    return value;
}

PropagatorGrid* Propagator::getGrid() const {
    return grid;
}

double Propagator::zeroPotential(double x) {
    return 0.0; //0.5 * x * x;
}

double Propagator::harmonicPotential(double x) {
    return 0.5 * x * x;
}

void Propagator::setPotential(double (*v)(double)) {
    potential = v;
}




