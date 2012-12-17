#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

class PropagatorGrid;

class Propagator {
public:
    Propagator();
    virtual ~Propagator();

    double evaluate();

    void setupGrid();
    void initializeGrid(int index0);
    void propagate();
    double readValue(int index) const;

    PropagatorGrid* getGrid() const;
private:
    PropagatorGrid* grid;

    static double potential(double x);

};

#endif
