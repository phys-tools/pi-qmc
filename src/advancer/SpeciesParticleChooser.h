#ifndef __SpeciesParticleChooser_h_
#define __SpeciesParticleChooser_h_

#include "ParticleChooser.h"
class Species;
class SimulationInfo;

/// Class for choosing particles at random from a given species.
class SpeciesParticleChooser : public ParticleChooser {
public:
  SpeciesParticleChooser(const Species&, const int nmoving);
  virtual ~SpeciesParticleChooser();

  virtual void chooseParticles();
protected:
  const int npart;
  const int ifirst;
};
#endif
