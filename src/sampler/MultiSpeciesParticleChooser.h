/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifndef __MultiSpeciesParticleChooser_h_
#define __MultiSpeciesParticleChooser_h_

#include "ParticleChooser.h"
class Species;
class SimulationInfo;

/// Class for choosing particles at random from a more than one species.

class MultiSpeciesParticleChooser : public ParticleChooser {
  typedef blitz::Array<int,1> IArray;
public:
  /// Construct by giving Species and number of moving particles.
  MultiSpeciesParticleChooser( const Species * speciesList, const int nspecies, const int nmoving);
  /// Virtual destructor.
  virtual ~MultiSpeciesParticleChooser();
  /// Choose a new set of particles.
  virtual void chooseParticles();

 protected:
  IArray speciesContainer;
  int count;
};
#endif
