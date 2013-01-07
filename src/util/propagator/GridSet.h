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

    void setupGrid();
    void initializeGrid();

    bool isConverged() const;

    double readValue0() const;
    PropagatorGrid* getGrid(int index);
private:
    typedef std::vector<PropagatorGrid*> GridArray;
    typedef std::vector<PropagatorGrid*>::iterator GridIterator;
    GridArray grid;
    GridParameters* parameters;
};

#endif
