#include "SimpleParticleChooser.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "util/RandomNumGenerator.h"
#include <iostream>

SimpleParticleChooser::SimpleParticleChooser(const int npart, const int nmoving) :
        ParticleChooser(nmoving), npart(npart) {
}

SimpleParticleChooser::~SimpleParticleChooser() {
}

void SimpleParticleChooser::chooseParticles() {
    for (int imoving = 0; imoving < index.size();) {
        int i = (int) (npart * RandomNumGenerator::getRand());
        if (i == npart)
            continue;
        for (int jmoving = 0; jmoving <= imoving; ++jmoving) {
            if (jmoving == imoving)
                index(imoving++) = i;
            if (index(jmoving) == i)
                break;
        }
    }
}
