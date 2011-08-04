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
#ifndef __ParticleChooser_h_
#define __ParticleChooser_h_

#include <cstdlib>
#include <blitz/array.h>
#include <string>
/// Base class for algorithms for choosing particles to move.
/// @version $Revision$
/// @author John Shumway
class ParticleChooser {
public:
  typedef blitz::Array<int,1> IArray;
  /// Construct by giving the number of particles to choose.
  ParticleChooser(const int npart);
  /// Virtual destructor.
  virtual ~ParticleChooser() {}
  /// Choose a new set of particles.
  virtual void chooseParticles()=0;
  /// Get the index of the ith moving particle.
  int operator[](const int i) const {return index(i);}
  /// Get name.
  std::string& getName() {return name;}
protected:
  /// The indicies of the particles to move.
  IArray index;
  /// A string identifier.
  std::string name;
};
#endif
