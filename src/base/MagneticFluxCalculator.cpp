#include "MagneticFluxCalculator.h"

#include "SimulationInfo.h"
#include "Paths.h"
#include "Charges.h"

MagneticFluxCalculator::MagneticFluxCalculator(Charges* charges)
:   charges(charges),
    flux(0.0) {
}

MagneticFluxCalculator::~MagneticFluxCalculator() {
    delete charges;
}

double MagneticFluxCalculator::calculate(Paths* paths) {
    paths->sumOverLinks(*this);
    return flux;
}

void MagneticFluxCalculator::initCalc(const int nslice, const int firstSlice) {
    flux = 0.0;
}

void MagneticFluxCalculator::handleLink(const Vec& start, const Vec& end,
        int ipart, int islice, const Paths& paths) {
    double startCrossEnd = start[0] * end[1] - start[1] * end[0];
    flux += charges->getValue(ipart) * 0.5 * startCrossEnd;
}

void MagneticFluxCalculator::endCalc(int nslice) {
}

