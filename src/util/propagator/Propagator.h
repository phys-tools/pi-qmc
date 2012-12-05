#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

class PropagatorGrid;

class Propagator {
public:
    Propagator();
    virtual ~Propagator();

    double evaluate();

    void setupGrid();
    void initializeGrid();
    void propagate();
    double readValue() const;

    PropagatorGrid* getGrid() const;
private:
    PropagatorGrid* grid;
};

#endif
