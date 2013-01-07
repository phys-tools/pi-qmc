#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

class PropagatorGrid;
class GridSet;

class Propagator {
public:
    Propagator(double mass, double tau, double x0);
    virtual ~Propagator();

    double evaluate();

    void propagate();

    double getGridSpacing() const;
    void setPotential(double (*v)(double));

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
};

#endif
