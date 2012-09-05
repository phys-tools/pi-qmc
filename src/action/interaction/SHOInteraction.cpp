#include "config.h"
#include "SHOInteraction.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"

SHOInteraction::SHOInteraction(const SimulationInfo& simInfo, double omega,
        const Species *species1, const Species *species2, int maxlevel) :
        omega(omega), mu(1.0 / (1.0 / species1->mass + 1.0 / species2->mass)), deltaTau(
                simInfo.getTau()), coshwt(maxlevel+1), sinhwt(maxlevel+1) {
    for (int i = 0; i <= maxlevel; ++i) {
        double deltaTauEff = deltaTau * (1 << i);
        coshwt(i) = cosh(omega * deltaTauEff);
        sinhwt(i) = sinh(omega * deltaTauEff);
    }
}

SHOInteraction::~SHOInteraction() {
}

double SHOInteraction::getActionDifference(
        const SectionSamplerInterface& sampler, int level) {
    double difference = 0.0;

    const Beads<NDIM> &movingBeads = sampler.getMovingBeads();
    const Beads<NDIM> &oldBeads = sampler.getSectionBeads();
    int nslice = movingBeads.getNSlice();

    Vec prevMovingDelta = movingBeads(0, 0) - movingBeads(1, 0);
    Vec prevOldDelta = oldBeads(0, 0) - oldBeads(1, 0);
    int step = (1 << level);
    for (int slice = step; slice < nslice; slice += step) {
        Vec nextMovingDelta = movingBeads(0, slice) - movingBeads(1, slice);
        Vec nextOldDelta = oldBeads(0, slice) - oldBeads(1, slice);

        difference += calculateAction(nextMovingDelta, prevMovingDelta, level);
        difference -= calculateAction(nextOldDelta, nextOldDelta, level);
        prevMovingDelta = nextMovingDelta;
        prevOldDelta = nextOldDelta;
    }
    return difference;
}

double SHOInteraction::getActionDifference(const Paths& paths,
        const VArray& displacement, int nmoving, const IArray& movingIndex,
        int iFirstSlice, int iLastSlice) {
    return 0.0;
}

double SHOInteraction::getTotalAction(const Paths& paths,
        const int level) const {
    return 0.0;
}

void SHOInteraction::getBeadAction(const Paths& paths, const int ipart,
        const int islice, double& u, double& utau, double& ulambda, Vec& fm,
        Vec& fp) const {
    if (ipart != 0) return;
    Vec delta1 = paths(0, islice, -1) - paths(1, islice, -1);
    Vec delta2 = paths(0, islice) - paths(1, islice);
    u = calculateAction(delta1, delta2, 0);
    utau = calculateTauDerivativeOfAction(delta1, delta2);
}

double SHOInteraction::calculateAction(Vec delta1, Vec delta2, int level) const {
    double dist1 = sqrt(dot(delta1, delta1));
    double dist2 = sqrt(dot(delta2, delta2));
    double deltaTauEff = deltaTau * (1 << level);
    double value = 0.5 * NDIM * log(sinhwt(level) / (omega * deltaTauEff));
    value += 0.5 * mu * omega
            * ((dist1 * dist1 + dist2 * dist2) * coshwt(level)
                    - 2.0 * dot(delta1, delta2)) / sinhwt(level);
    value -= 0.5 * mu * (dist1 - dist2) * (dist1 - dist2) / deltaTauEff;
    return value;
}

double SHOInteraction::calculateTauDerivativeOfAction(Vec delta1,
        Vec delta2) const {
    double dist1 = sqrt(dot(delta1, delta1));
    double dist2 = sqrt(dot(delta2, delta2));
    double value = 0.5 * NDIM * omega
            * (coshwt(0) / sinhwt(0) - 1.0 / (omega * deltaTau));
    value += 0.5 * mu * omega * omega * (dist1 * dist1 + dist2 * dist2);
    value -= 0.5 * mu * omega * coshwt(0)  / (sinhwt(0) * sinhwt(0))
            * ((dist1 * dist1 + dist2 * dist2) * coshwt(0)
                    - 2.0 * dot(delta1, delta2));
    value += 0.5 * mu * (dist1 - dist2) * (dist1 - dist2)
            / (deltaTau * deltaTau);
    return value;
}


