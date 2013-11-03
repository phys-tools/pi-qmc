#ifndef __TwoPairChooser_h_
#define __TwoPairChooser_h_

#include "ParticleChooser.h"
#include "PermutationChooser.h"
class Species;

class TwoPairChooser: public ParticleChooser, public PermutationChooser {
public:
    TwoPairChooser(const Species&, const Species&);
    virtual ~TwoPairChooser();

    virtual void chooseParticles();
private:
    const int npart1, npart2, ifirst1, ifirst2;
};
#endif
