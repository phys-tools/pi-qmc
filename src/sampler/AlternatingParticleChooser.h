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
#ifndef __AlternatingParticleChooser_h_
#define __AlternatingParticleChooser_h_

#include "ParticleChooser.h"
class Species;
class SimulationInfo;

/// Class for choosing particles from different species in an alternating way.

class AlternatingParticleChooser : public ParticleChooser {
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
