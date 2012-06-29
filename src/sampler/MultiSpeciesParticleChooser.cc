#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MultiSpeciesParticleChooser.h"
#include "util/RandomNumGenerator.h"
#include "SimulationInfo.h"
#include <iostream>

MultiSpeciesParticleChooser::MultiSpeciesParticleChooser(
        const Species *speciesList, const int nspecies, const int nmoving) :
        ParticleChooser(nmoving) {

    count = 0;
    for (int ispec = 0; ispec < nspecies; ispec++) {
        count += speciesList[ispec].count;
    }

    speciesContainer.resize(count);
    int k = 0;
    for (int ispec = 0; ispec < nspecies; ispec++) {
        for (int i = 0; i < speciesList[ispec].count; i++) {
            speciesContainer(k++) = speciesList[ispec].ifirst + i;
        }
    }
}

MultiSpeciesParticleChooser::~MultiSpeciesParticleChooser() {
}

void MultiSpeciesParticleChooser::chooseParticles() {
    for (int imoving = 0; imoving < index.size(); imoving++) {
        int i = (int) (count * RandomNumGenerator::getRand());
        index(imoving) = speciesContainer(i);
    }
}
