#ifndef POTENTIALGRID_H_
#define POTENTIALGRID_H_

class PotentialGrid {
public:
    struct functor {virtual double operator()(double) {return 0.0;}};

    PotentialGrid(int size, double deltaX, double x0, double (*v)(double),
                  double deltaTau);
    PotentialGrid(int size, double deltaX, double x0,
                  functor& v,
                  double deltaTau);
    virtual ~PotentialGrid();

    double operator()(int index) const;
private:
    double* value;
    int size;
};

#endif
