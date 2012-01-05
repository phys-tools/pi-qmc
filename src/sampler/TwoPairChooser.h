// $Id: TwoPairChooser.h 338 2010-11-30 18:56:16Z john.shumwayjr $
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
#ifndef __TwoPairChooser_h_
#define __TwoPairChooser_h_

#include "ParticleChooser.h"
#include "PermutationChooser.h"
class Species;
/// Class for choosing  particles at random.
/// @version $Revision: 338 $
/// @author John Shumway
class TwoPairChooser : public ParticleChooser, public PermutationChooser {
public:
  /// Construct by giving the two species.
  TwoPairChooser(const Species&, const Species&);
  /// Virtual destructor.
  virtual ~TwoPairChooser();
  /// Choose a new set of particles.
  virtual void chooseParticles();
private:
  /// The number of particles.
  const int npart1,npart2,ifirst1,ifirst2;
};
#endif
