#ifndef __MultiSpeciesParticleChooser_h_
#define __MultiSpeciesParticleChooser_h_

#include "ParticleChooser.h"
class Species;
class SimulationInfo;

/// Class for choosing particles at random from a more than one species.
class MultiSpeciesParticleChooser: public ParticleChooser {
    typedef blitz::Array<int, 1> IArray;
public:
    MultiSpeciesParticleChooser(const Species * speciesList, const int nspecies,
            const int nmoving);

    virtual ~MultiSpeciesParticleChooser();

    virtual void chooseParticles();

protected:
    IArray speciesContainer;
    int count;
};
#endif
