// $Id$
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
#ifndef __SpeciesParticleChooser_h_
#define __SpeciesParticleChooser_h_

#include "ParticleChooser.h"
class Species;
class SimulationInfo;

/// Class for choosing particles at random from a given species.
/// @version $Revision$
/// @author John Shumway  
class SpeciesParticleChooser : public ParticleChooser {
public:
  /// Construct by giving Species and number of moving particles.
  SpeciesParticleChooser(const Species&, const int nmoving);
  /// Virtual destructor.
  virtual ~SpeciesParticleChooser();
  /// Choose a new set of particles.
  virtual void chooseParticles();
protected:
  /// The number of particles in this species.
  const int npart;
  /// The index of the first particle in the species.
  const int ifirst;
};
#endif
