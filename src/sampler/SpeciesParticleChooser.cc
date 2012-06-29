#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "SpeciesParticleChooser.h"
#include "util/RandomNumGenerator.h"
#include "SimulationInfo.h"
#include <iostream>

SpeciesParticleChooser::SpeciesParticleChooser(const Species& species,
        const int nmoving) :
        ParticleChooser(nmoving), npart(species.count), ifirst(species.ifirst) {
    name = species.name;
}

SpeciesParticleChooser::~SpeciesParticleChooser() {
}

void SpeciesParticleChooser::chooseParticles() {
    for (int imoving = 0; imoving < index.size();) {
        int i = (int) (npart * RandomNumGenerator::getRand() + ifirst);
        if (i == npart + ifirst)
            i = npart + ifirst - 1;
        for (int jmoving = 0; jmoving <= imoving; ++jmoving) {
            if (jmoving == imoving)
                index(imoving++) = i;
            if (index(jmoving) == i)
                break;
        }
    }
}
