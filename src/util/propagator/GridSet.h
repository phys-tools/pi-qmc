#ifndef GRIDSET_H_
#define GRIDSET_H_

#include <vector>
class PropagatorGrid;
class GridParameters;

class GridSet {
public:
    GridSet();
    virtual ~GridSet();

    double getDeltaX() const;
    void setGridParameters(GridParameters*);

    void extrapolateValue0();
    bool isConverged() const;
    double readValue0() const;

    PropagatorGrid* makeNewGrid(double deltaTau);
private:
    typedef std::vector<PropagatorGrid*> GridArray;
    typedef std::vector<double> Array;
    GridArray gridArray;
    GridParameters* parameters;
    Array deltaTau;
    PropagatorGrid* allocateGrid();
    int gridCount;
    double solution;
    double lastDelta;
    static const double TOLERANCE;
};

#endif
