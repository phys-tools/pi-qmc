#ifndef PROPAGATOR3D_H_
#define PROPAGATOR3D_H_

#include "Propagator1D.h"

class Propagator3D {
public:
  Propagator3D(double mass, double tau, double x0);
  virtual ~Propagator3D();

  double evaluate();

  void propagate();

  double getGridSpacing() const;
  void setPotential(double (*v)(double));
  void setPotential(PotentialGrid::functor*);

private:
  Propagator1D prop_1D;
};

#endif
