#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "PairChooser.h"
#include "WalkingChooser.h"
#include "util/Permutation.h"
#include <iostream>

PairChooser::PairChooser(const int npart, const Species &s1, const Species &s2,
        const int nlevel, const SimulationInfo &simInfo) :
        ParticleChooser(2 * npart), PermutationChooser(2 * npart), npart(npart), chooser1(
                *new WalkingChooser(npart, s1, nlevel, simInfo)), chooser2(
                *new WalkingChooser(npart, s2, nlevel, simInfo)) {
}

PairChooser::~PairChooser() {
    delete &chooser1;
    delete &chooser2;
}

void PairChooser::chooseParticles() {
}

bool PairChooser::choosePermutation() {
    bool newPerm = chooser1.choosePermutation() && chooser2.choosePermutation();
    if (newPerm) {
        for (int i = 0; i < npart; ++i) {
            index(i) = chooser1[i];
            index(i + npart) = chooser2[i];
            (*permutation)[i] = chooser1.getPermutation()[i];
            (*permutation)[i + npart] = chooser2.getPermutation()[i] + npart;
        }
    }
    return newPerm;
}

double PairChooser::getLnTranProb() const {
    return chooser1.getLnTranProb() + chooser2.getLnTranProb();
}

void PairChooser::init() {
    chooser1.init();
    chooser2.init();
}

void PairChooser::setMLSampler(const MultiLevelSampler *sampler) {
    chooser1.setMLSampler(sampler);
    chooser2.setMLSampler(sampler);
}
