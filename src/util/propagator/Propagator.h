#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

class PropagatorGrid;

class Propagator {
public:
    Propagator(double tau);
    virtual ~Propagator();

    double evaluate();

    void setupGrid();
    void initializeGrid(int index0);
    void propagate();
    double readValue(int index) const;

    PropagatorGrid* getGrid() const;
    void setPotential(double (*v)(double));

    static double zeroPotential(double x);
    static double harmonicPotential(double x);
private:
    PropagatorGrid* grid;
    double tau;
    double (*potential)(double);


};

#endif
