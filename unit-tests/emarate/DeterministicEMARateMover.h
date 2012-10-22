#ifndef DETERMINISTICEMARATEMOVER_H_
#define DETERMINISTICEMARATEMOVER_H_

#include "emarate/EMARateMover.h"

class DeterministicEMARateMover : public EMARateMover {
public:
    DeterministicEMARateMover(double tau,
            const Species *species1, const Species *species2,
            int maxlevel, double C)
    : EMARateMover(tau, species1, species2, maxlevel, C),
      nextRandomNumber(0.5) {
    }

    virtual double getRandomNumber() const {
        return nextRandomNumber;
    }

    virtual void makeGaussianRandomNumbers(blitz::Array<Vec,1> &gaussRand) {
        gaussRand = nextGaussianRandomNumber;
    }

    double nextRandomNumber;
    Vec nextGaussianRandomNumber;
};
#endif
