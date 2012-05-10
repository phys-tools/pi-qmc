#ifndef DETERMINISTICEMARATEMOVER_H_
#define DETERMINISTICEMARATEMOVER_H_

#include "emarate/EMARateMover.h"

class DeterministicEMARateMover : public EMARateMover {
public:
    DeterministicEMARateMover(double tau, double mass1, double mass2,
            int maxlevel, double C)
    : EMARateMover(tau, mass1, mass2, maxlevel, C),
      nextRandomNumber(0.5) {
    }

    virtual double getRandomNumber() const {
        return nextRandomNumber;
    }

    double nextRandomNumber;
};
#endif
