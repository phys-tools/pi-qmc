#ifndef LENNARDJONESPOTENTIAL_H_
#define LENNARDJONESPOTENTIAL_H_

#include "PairPotential.h"

/// Lennard-Jones 6-12.
class LennardJonesPotential: public PairPotential {
public:
    LennardJonesPotential(double epsilon, double sigma) :
            epsilon(epsilon), sigma(sigma) {
    }
    virtual ~LennardJonesPotential() {}
    virtual double operator()(double r) const {
        double x = sigma / r;
        x *= x * x;
        x *= x;
        return 4 * epsilon * x * (x - 1);
    }
    double epsilon, sigma;
};

#endif
