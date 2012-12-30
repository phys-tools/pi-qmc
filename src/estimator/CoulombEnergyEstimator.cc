#include "config.h"
#include "CoulombEnergyEstimator.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "base/Paths.h"
#include "stats/ScalarAccumulator.h"
#include "util/SuperCell.h"
#include <blitz/tinyvec-et.h>

CoulombEnergyEstimator::CoulombEnergyEstimator(const SimulationInfo& simInfo,
        const double epsilon, const std::string& unitName,
        double scale, double shift, ScalarAccumulator *accumulator)
    :   ScalarEstimator("coulomb_energy", "scalar-energy/coulomb-energy",
                    unitName, scale, shift),
        epsilon(epsilon),
        q(simInfo.getNPart()) {
    this->accumulator = accumulator;
    for (int i = 0; i < q.size(); ++i) {
        q(i) = simInfo.getPartSpecies(i).charge;
    }
}

CoulombEnergyEstimator::~CoulombEnergyEstimator() {
}

void CoulombEnergyEstimator::initCalc(const int nslice, const int firstSlice) {
    accumulator->clearValue();
}

void CoulombEnergyEstimator::handleLink(const Vec& start, const Vec& end,
        const int ipart, const int islice, const Paths& paths) {
    const SuperCell& cell = paths.getSuperCell();
    double energy = 0.0;
    for (int jpart = 0; jpart < ipart; ++jpart) {
        Vec delta = end - paths(jpart, islice);
        cell.pbc(delta);
        energy += q(ipart) * q(jpart) / (epsilon * sqrt(dot(delta, delta)));
    }
    accumulator->addToValue(energy);
}

void CoulombEnergyEstimator::endCalc(const int lnslice) {
    accumulator->storeValue(lnslice);
}

double CoulombEnergyEstimator::calcValue() {
    return 0.0;
}

void CoulombEnergyEstimator::reset() {
    accumulator->reset();
}

void CoulombEnergyEstimator::evaluate(const Paths& paths) {
    paths.sumOverLinks(*this);
}

