#ifndef INVERSECOSH2POTENTIAL_H_
#define INVERSECOSH2POTENTIAL_H_
#include "PairPotential.h"
#include <cmath>

/// Empirical short range potential commonly used for trapped atoms.
class InverseCosh2Potential : public PairPotential {
public:
    InverseCosh2Potential(double v0, double kappa) :
            v0(v0), kappa(kappa) {
    }
    virtual ~InverseCosh2Potential() {}
    virtual double operator()(double r) const {
        double coshKappaR = cosh(kappa * r);
        return v0 / (coshKappaR * coshKappaR);
    }
    double v0, kappa;
};
#endif
