#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "TwoPairChooser.h"
#include "base/Species.h"
#include "util/RandomNumGenerator.h"
#include "util/Permutation.h"
#include <iostream>

TwoPairChooser::TwoPairChooser(const Species &s1, const Species &s2) :
        ParticleChooser(4), PermutationChooser(4), npart1(s1.count), npart2(
                s2.count), ifirst1(s1.ifirst), ifirst2(s2.ifirst) {
    (*this->permutation)[0] = 1;
    (*this->permutation)[1] = 0;
    (*this->permutation)[2] = 3;
    (*this->permutation)[3] = 2;
}

TwoPairChooser::~TwoPairChooser() {
}

void TwoPairChooser::chooseParticles() {
    int i = 0, j = 0, k = 0, l = 0;
    do {
        i = (int) (npart1 * RandomNumGenerator::getRand()) + ifirst1;
    } while (i == npart1 + ifirst1);
    do {
        j = (int) (npart1 * RandomNumGenerator::getRand()) + ifirst1;
    } while (j == i || j == npart1 + ifirst1);
    do {
        k = (int) (npart2 * RandomNumGenerator::getRand()) + ifirst2;
    } while (k == i || k == j || k == npart2 + ifirst2);
    do {
        l = (int) (npart2 * RandomNumGenerator::getRand()) + ifirst2;
    } while (l == k || l == i || l == j || l == npart2 + ifirst2);
    index(0) = i;
    index(1) = j;
    index(2) = k;
    index(3) = l;
}
