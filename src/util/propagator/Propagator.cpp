#include "Propagator.h"
#include "PropagatorGrid.h"
#include "GridParameters.h"
#include "GridSet.h"

Propagator::Propagator(double mass, double tau, double x0)
    :   tau(tau),
        mass(mass),
        x0(x0),
        deltaX(0.005) {
    gridSet = new GridSet();
    gridSet->setGridParameters(new GridParameters(mass, tau, x0, deltaX));
    potential = zeroPotential;
}

Propagator::~Propagator() {
    delete gridSet;
}

double Propagator::evaluate() {
    propagate();
    double value = gridSet->readValue0();
    return value;
}

void Propagator::propagate() {
    int istep = 0;
    propagate(istep++);
    propagate(istep++);
    while (! gridSet->isConverged()) {
        propagate(istep++);
    }
}

void Propagator::propagate(int step) {
    int stepCount = 2 * (step + 1);
    double deltaTau = tau / stepCount;
    PropagatorGrid* grid = gridSet->makeNewGrid(deltaTau);
    propagate(grid, deltaTau, stepCount);
}

void Propagator::propagate(PropagatorGrid* grid, double deltaTau,
        int stepCount) {
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
    return gridSet->getDeltaX();
}

