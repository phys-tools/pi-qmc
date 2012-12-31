#ifndef GRIDPARAMETERS_H_
#define GRIDPARAMETERS_H_

class GridParameters {
public:
    GridParameters();
    virtual ~GridParameters();

    int getGridCount() const;
    double getDeltaX() const;
    double getXMin() const;
    int getIndex0() const;
private:
    int gridCount;
    double deltaX;
    double xmin;
    int index0;
};

#endif
