#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "EMARateWeight.h"
#include "action/coulomb/CoulombLinkAction.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "base/Paths.h"
#include "util/SuperCell.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

EMARateWeight::EMARateWeight(const SimulationInfo& simInfo,
        const Species* species1, const Species* species2, double C)
:   dtau(simInfo.getTau()),
    C(C),
    cell(*simInfo.getSuperCell()),
    lastSlice(simInfo.getNSlice() - 1),
    sum(0),
    norm(0),
    hasCoulomb(false),
    coulomb(0),
    species1(species1),
    species2(species2),
    index1(species1->ifirst),
    index2(species2->ifirst) {
    if (species1->anMass) {
        mass1 = *(species1->anMass);
    } else {
        mass1 = species1->mass;
    }
    if (species2->anMass) {
        mass2 = *(species2->anMass);
    } else {
        mass2 = species2->mass;
    }
}

EMARateWeight::~EMARateWeight() {
    delete coulomb;
}

void EMARateWeight::includeCoulombContribution(double epsilon, int norder) {
    hasCoulomb = true;
    double q1q2 = -1.0;
    double mu = 1.0 / (1.0 / (species1->mass) + 1.0 / (species2->mass));
    coulomb = new CoulombLinkAction(q1q2, epsilon, mu, dtau, norder);
}


void EMARateWeight::evaluateElectronBeforeRecombination(const Vec &start,
        const Vec &end, const Paths &paths) {
    Vec delta = paths(index1, 0) - start;
    cell.pbc(delta);
    for (int idim = 0; idim < NDIM; ++idim) {
        actionDifference +=
                0.5 * mass2(idim) * delta(idim)  * delta(idim) / dtau;
    }

    delta = end - start;
    cell.pbc(delta);
    for (int idim = 0; idim < NDIM; ++idim) {
        actionDifference -=
                0.5 * mass2(idim) * delta(idim) * delta(idim) / dtau;
    }
    if (hasCoulomb) {
        Vec zero = 0.0;
        Vec deltaPrev = paths(index1, 0, -1) - start;
        actionDifference += coulomb->getValue(deltaPrev, zero);

        delta = paths(index1, 0) - end;
        actionDifference -= coulomb->getValue(deltaPrev, delta);
    }
}

int EMARateWeight::getPartitionCount() const {
    return 2;
}

void EMARateWeight::evaluate(Paths* paths) {
    actionDifference = 0.0;
    Vec start = (*paths)(index2, 0, -1);
    Vec end = (*paths)(index2, 0);
    evaluateElectronBeforeRecombination(start, end, *paths);
    start = (*paths)(index1, 1, -1);
    end = (*paths)(index1, 1);
    evaluateHoleAfterRecombination(start, end, *paths);
    if (C > 0.0 && (actionDifference - log(C)) > -40.0) {
        diagonalProbability = 1.0 / (1.0 + C * exp(-actionDifference));
        recombiningProbability = 1.0 - diagonalProbability;
    } else {
        diagonalProbability = 0.0;
        recombiningProbability = 1.0;
    }

}

double EMARateWeight::getValue(int i) const {
    if (i==0) {
        return diagonalProbability;
    } else {
        return recombiningProbability;
    }
}

void EMARateWeight::evaluateHoleAfterRecombination(const Vec &start,
        const Vec &end, const Paths &paths) {
    Vec delta = end - paths(index2, 1, -1);
    cell.pbc(delta);
    for (int idim = 0; idim < NDIM; ++idim) {
        actionDifference +=
                0.5 * mass1(idim) * delta(idim) * delta(idim) / dtau;
    }
    delta = end - start;
    cell.pbc(delta);
    for (int idim = 0; idim < NDIM; ++idim) {
        actionDifference -=
                0.5 * mass1(idim) * delta(idim) * delta(idim) / dtau;
    }
    if (hasCoulomb) {
        Vec zero = 0.0;
        Vec deltaNext = paths(index2, 1) - end;
        actionDifference += coulomb->getValue(deltaNext, zero);

        delta = paths(index2, 1, -1) - start;
        actionDifference -= coulomb->getValue(deltaNext, delta);
    }
}
