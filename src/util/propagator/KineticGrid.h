#ifndef KINETICGRID_H_
#define KINETICGRID_H_

class KineticGrid {
public:
    KineticGrid(int gridCount, double deltaK, double mass, double deltaTau);
    virtual ~KineticGrid();

    double operator() (int index) const;
private:
    int gridCount;
    double *value;
};

#endif
