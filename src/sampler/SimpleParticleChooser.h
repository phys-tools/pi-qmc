#ifndef __SimpleParticleChooser_h_
#define __SimpleParticleChooser_h_

#include "ParticleChooser.h"

/// Class for choosing  particles at random.
class SimpleParticleChooser: public ParticleChooser {
public:
    SimpleParticleChooser(const int npart, const int nmoving);
    virtual ~SimpleParticleChooser();

    virtual void chooseParticles();
private:
    const int npart;
};
#endif
