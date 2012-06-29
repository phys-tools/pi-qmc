#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "ParticleChooser.h"
ParticleChooser::ParticleChooser(const int npart) :
        index(npart) {
    for (int i = 0; i < npart; ++i)
        index(i) = i;
}
