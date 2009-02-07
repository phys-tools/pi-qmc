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
#ifndef __SimpleParticleChooser_h_
#define __SimpleParticleChooser_h_

#include "ParticleChooser.h"
/// Class for choosing  particles at random.
/// @version $Revision$
/// @author John Shumway
class SimpleParticleChooser : public ParticleChooser {
public:
  /// Construct by giving the number of particles and number of moving particles.
  SimpleParticleChooser(const int npart, const int nmoving);
  /// Virtual destructor.
  virtual ~SimpleParticleChooser();
  /// Choose a new set of particles.
  virtual void chooseParticles();
private:
  /// The number of particles.
  const int npart;
};
#endif
