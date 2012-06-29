#ifndef __AlternatingParticleChooser_h_
#define __AlternatingParticleChooser_h_

#include "ParticleChooser.h"
class Species;
class SimulationInfo;

/// Class for choosing particles from different species in an alternating way.
class AlternatingParticleChooser: public ParticleChooser {
public:
    /// Construct by giving 2 Species.
    AlternatingParticleChooser(const Species&, const Species&, const int);
    /// Virtual destructor.
    virtual ~AlternatingParticleChooser();
    /// Choose a pair of particles.
    virtual void chooseParticles();
protected:
    /// The number of particles in each species.
    const int npart1, npart2;
    /// The index of the first particle in each species.
    const int ifirst1, ifirst2;
};
#endif
