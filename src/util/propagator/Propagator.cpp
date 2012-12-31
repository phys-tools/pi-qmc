#include "Propagator.h"
#include "PropagatorGrid.h"
#include "GridParameters.h"

Propagator::Propagator(double tau, double mass)
    :   grid(0),
        tau(tau),
        mass(mass) {
    potential = zeroPotential;
    gridParameters = new GridParameters();
}

Propagator::~Propagator() {
    delete grid;
    delete gridParameters;
}

double Propagator::evaluate() {
    setupGrid();
    int index0 = gridParameters->getIndex0();
    initializeGrid(index0);
    propagate();
    double value = readValue(index0);
    return value;
}

void Propagator::setupGrid() {
    int gridCount = gridParameters->getGridCount();
    double deltaX = gridParameters->getDeltaX();
    double xmin = gridParameters->getXMin();
    grid = new PropagatorGrid(gridCount, deltaX, xmin);
}

void Propagator::initializeGrid(int index0) {
    grid->initialize(index0);
}

void Propagator::propagate() {
    int stepCount = 100;
    double deltaTau = tau / stepCount;
    mass = 1.0;
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
    return 0.0;
}

double Propagator::harmonicPotential(double x) {
    return 0.5 * x * x;
}

void Propagator::setPotential(double (*v)(double)) {
    potential = v;
}

double Propagator::getGridSpacing() const {
    grid->getDeltaX();
}





