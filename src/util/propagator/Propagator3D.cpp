#include "Propagator3D.h"

Propagator3D::Propagator3D(double mass, double tau, double x0)
  : prop_1D(mass, tau, x0)
{
}

Propagator3D::~Propagator3D() {
}

double Propagator3D::evaluate() {
  return prop_1D.evaluate();
}

void Propagator3D::propagate() {
  prop_1D.propagate();
}

double Propagator3D::getGridSpacing() const {
  return prop_1D.getGridSpacing();
}

void Propagator3D::setPotential(double (*v)(double)) {
  prop_1D.setPotential(v);
}

void Propagator3D::setPotential(PotentialGrid::functor *v) {
  prop_1D.setPotential(v);
}
