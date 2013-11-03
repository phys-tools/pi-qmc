#ifndef POTENTIALGRID_H_
#define POTENTIALGRID_H_

class PotentialGrid {
public:
    PotentialGrid(int size, double deltaX, double x0, double (*v)(double),
            double deltaTau);
    virtual ~PotentialGrid();

    double operator()(int index) const;
private:
    double* value;
    int size;
};

#endif
