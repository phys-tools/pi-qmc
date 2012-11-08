#ifndef AZIZPOTENTIAL_H_
#define AZIZPOTENTIAL_H_

#include "PairPotential.h"

/// Aziz potential (HFD-B2).
/// R. A. Aziz et al., Mol. Phys. 77, 321 (1992).
class AzizPotential: public PairPotential {
public:
    AzizPotential() {}
    virtual ~AzizPotential() {}
    virtual double operator()(double r) const {
        double x = r / r0;
        double v = (x < D) ? exp(-(D / x - 1) * (D / x - 1)) : 1.;
        double x2inv = 1. / (x * x);
        v *= -x2inv * x2inv * x2inv * (C6 + x2inv * (C8 + x2inv * C10));
        v += A * exp(-x * (alpha - beta * x));
        v *= eps;
        return v;
    }
    static const double A;
    static const double alpha;
    static const double beta;
    static const double C6;
    static const double C8;
    static const double C10;
    static const double D;
    static const double eps;
    static const double r0;
};

#endif
