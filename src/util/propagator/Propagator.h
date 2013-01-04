#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

class PropagatorGrid;
class GridParameters;

class Propagator {
public:
    Propagator(double mass, double tau, double x0);
    virtual ~Propagator();

    double evaluate();

    void setupGrid();
    void initializeGrid(int index0);
    void propagate();
    double readValue(int index) const;
    double getGridSpacing() const;

    PropagatorGrid* getGrid() const;
    void setPotential(double (*v)(double));

    static double zeroPotential(double x);
    static double harmonicPotential(double x);
private:
    PropagatorGrid* grid;
    GridParameters* gridParameters;
    double tau;
    double mass;
    double x0;
    double deltaX;
    double (*potential)(double);
};

#endif
