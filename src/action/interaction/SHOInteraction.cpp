#include "config.h"
#include "SHOInteraction.h"

SHOInteraction::SHOInteraction(const SimulationInfo& simInfo, double omega) {
}

SHOInteraction::~SHOInteraction() {
}

double SHOInteraction::getActionDifference(
        const SectionSamplerInterface& sampler, int level) {
    return 0.0;
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
