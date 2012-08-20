#include "config.h"
#include "SHOInteraction.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"

SHOInteraction::SHOInteraction(const SimulationInfo& simInfo, double omega,
        const Species *species1, const Species *species2) :
        omega(omega), mu(1.0 / (1.0 / species1->mass + 1.0 / species2->mass)), deltaTau(
                simInfo.getTau()), coshwt(cosh(omega * deltaTau)), sinhwt(
                sinh(omega * deltaTau)) {
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
    for (int islice = 1; islice < nslice; ++islice) {
        Vec nextMovingDelta = movingBeads(0, islice) - movingBeads(1, islice);
        Vec nextOldDelta = oldBeads(0, islice) - oldBeads(1, islice);

        difference += calculateAction(nextMovingDelta, prevMovingDelta);
        difference -= calculateAction(nextOldDelta, nextOldDelta);
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
}
double SHOInteraction::calculateAction(Vec delta1, Vec delta2) {
    double dist1 = sqrt(dot(delta1, delta1));
    double dist2 = sqrt(dot(delta2, delta2));
    double value = 0.5 * log(sinhwt / (omega * deltaTau));
    value += 0.5 * mu * omega
            * ((dist1 * dist1 + dist2 * dist2) * coshwt
                    - 2.0 * dot(delta1, delta2)) / sinhwt;
    value -= 0.5 * mu * (dist1 - dist2) * (dist1 - dist2) / deltaTau;
    return value;
}

