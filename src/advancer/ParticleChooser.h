#ifndef __ParticleChooser_h_
#define __ParticleChooser_h_

#include <cstdlib>
#include <blitz/array.h>
#include <string>

/// Base class for algorithms for choosing particles to move.
class ParticleChooser {
public:
    typedef blitz::Array<int, 1> IArray;

    ParticleChooser(const int npart);

    virtual ~ParticleChooser() {
    }

    virtual void chooseParticles()=0;

    int operator[](const int i) const {
        return index(i);
    }

    std::string& getName() {
        return name;
    }
protected:
    IArray index;
    std::string name;
};
#endif
