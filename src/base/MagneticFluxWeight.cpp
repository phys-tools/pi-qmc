#include "MagneticFluxWeight.h"
#include "MagneticFluxCalculator.h"
#include <cmath>

MagneticFluxWeight::MagneticFluxWeight(double bfield, int partitionCount,
        MagneticFluxCalculator* fluxCalculator)
:   deltaB(bfield / (partitionCount - 1)),
    deltaPhi(0.0),
    partitionCount(partitionCount),
    fluxCalculator(fluxCalculator) {
}

MagneticFluxWeight::~MagneticFluxWeight() {
    delete fluxCalculator;
}

int MagneticFluxWeight::getPartitionCount() const {
    return partitionCount;
}

void MagneticFluxWeight::evaluate(Paths* paths) {
    double flux = fluxCalculator->calculate(paths);
    deltaPhi = flux * deltaB;
}

double MagneticFluxWeight::getValue(int i) const {
    return cos(deltaPhi * i);
}

