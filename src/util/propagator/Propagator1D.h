#ifndef PROPAGATOR1D_H_
#define PROPAGATOR1D_H_

#include "PotentialGrid.h"

class PropagatorGrid;
class GridSet;

class Propagator1D {
public:
    Propagator1D(double mass, double tau, double x0);
    virtual ~Propagator1D();

    double evaluate();

    void propagate();

    double getGridSpacing() const;
    void setPotential(double (*v)(double));
    void setPotential(PotentialGrid::functor*);

    static double zeroPotential(double x);
    static double harmonicPotential(double x);

private:
    GridSet* gridSet;
    double tau;
    double mass;
    double x0;
    double deltaX;
    void propagate(int step);
    void propagate(PropagatorGrid* grid, double deltaTau, int stepCount);
    double (*potential)(double);
    PotentialGrid::functor* potential_functor;
};

#endif
