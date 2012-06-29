#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "AlternatingParticleChooser.h"
#include "util/RandomNumGenerator.h"
#include "SimulationInfo.h"
#include <iostream>

AlternatingParticleChooser::AlternatingParticleChooser(const Species& species1,
        const Species& species2, const int nmoving) :
        ParticleChooser(nmoving), npart1(species1.count), npart2(
                species2.count), ifirst1(species1.ifirst), ifirst2(
                species2.ifirst) {
}

AlternatingParticleChooser::~AlternatingParticleChooser() {
}

void AlternatingParticleChooser::chooseParticles() {
    int i = (int) (npart1 * RandomNumGenerator::getRand() + ifirst1);
    if (i == ifirst1 + npart1)
        i = ifirst1 + npart1 - 1;
    index(0) = i;
    i = (int) (npart2 * RandomNumGenerator::getRand() + ifirst2);
    if (i == ifirst2 + npart2)
        i = ifirst2 + npart2 - 1;
    index(1) = i;
    if (RandomNumGenerator::getRand() < 0.5) {
        int temp = index(0);
        index(0) = index(1);
        index(1) = temp;
    }
//  std::cout<<"index(0) = "<<index(0)<<", index(1) = "<<index(1)<<std::endl;
//  std::cin.ignore();
}

