#ifndef GRIDPARAMETERS_H_
#define GRIDPARAMETERS_H_

class GridParameters {
public:
    GridParameters();
    GridParameters(double mass, double tau, double x0, double deltaX);
    virtual ~GridParameters();

    int getGridCount() const;
    double getDeltaX() const;
    double getXMin() const;
    int getIndex0() const;

    static double calculateThermalWidth(double mass, double tau);
    static int powerTwoCeiling(double width);
private:
    int gridCount;
    double deltaX;
    double xmin;
    int index0;
};

#endif
