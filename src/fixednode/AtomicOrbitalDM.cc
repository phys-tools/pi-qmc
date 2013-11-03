#include "AtomicOrbitalDM.h"

AtomicOrbitalDM::AtomicOrbitalDM(
    int ifirst, int npart, int nfermion, double weight)
  : weight(weight), ifirst(ifirst), npart(npart),
    dist1(nfermion,npart), dist2(nfermion,npart) {
}
