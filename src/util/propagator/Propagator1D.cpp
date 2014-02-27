#include "Propagator1D.h"
#include "PropagatorGrid.h"
#include "GridParameters.h"
#include "GridSet.h"

Propagator1D::Propagator1D(double mass, double tau, double x0)
:   tau(tau),
    mass(mass),
    x0(x0),
    deltaX(0.005),
    potential_functor(0) {
  gridSet = new GridSet();
  gridSet->setGridParameters(new GridParameters(mass, tau, x0, deltaX));
  potential = zeroPotential;
}

Propagator1D::~Propagator1D() {
  delete gridSet;
}

double Propagator1D::evaluate() {
  propagate();
  double value = gridSet->readValue0();
  return value;
}

void Propagator1D::propagate() {
  int inc = 5;
  int istep = 0;
  propagate(istep += inc);
  propagate(istep += inc);
  gridSet->extrapolateValue0();
  while (! gridSet->isConverged()) {
    propagate(istep += inc);
    gridSet->extrapolateValue0();
  }
}

void Propagator1D::propagate(int step) {
  int stepCount = 2 * (step + 1);
  double deltaTau = tau / stepCount;
  PropagatorGrid* grid = gridSet->makeNewGrid(deltaTau);
  propagate(grid, deltaTau, stepCount);
}

void Propagator1D::propagate(PropagatorGrid* grid, double deltaTau,
    int stepCount) {
  grid->setupKineticPropagator(mass, deltaTau);
  if (potential_functor)
    grid->setupPotentialPropagator(*potential_functor, deltaTau);
  else
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


double Propagator1D::zeroPotential(double x) {
  return 0.0;
}

double Propagator1D::harmonicPotential(double x) {
  return 0.5 * x * x;
}

void Propagator1D::setPotential(double (*v)(double)) {
  potential = v;
}

void Propagator1D::setPotential(PotentialGrid::functor* v) {
  potential_functor = v;
}

double Propagator1D::getGridSpacing() const {
  return gridSet->getDeltaX();
}

