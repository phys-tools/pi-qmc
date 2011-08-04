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
#ifndef __Species_h_
#define __Species_h_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <iostream>
#include <cstdlib>
#include <blitz/tinyvec.h>

/** Class for storing species information.
 * @version $Revision$
 * @author John Shumway */
class Species {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Default constructor.
  Species();
  /// Constructor.
  Species(const std::string& name, const int count, const double mass,
          const double charge, const int twos, const bool isFermion);
  /// Copy constructor.
  Species(const Species& s);
  /// Copy operator.
  Species& operator=(const Species& s);
  /// Destructor.
  ~Species();
  /// The name of the species.
  std::string name;
  /// The number of particles.
  int count;
  /// The index of the first particle of this species.
  int ifirst;
  /// The mass of the species.
  double mass;
  /// The charge of the species.
  double charge;
  /// The spin of the species (in units of 1/2).
  int twos;
  /// Flag for whether the particles have fermi statistics.
  bool isFermion;
  /// Flag to indicate whether particle is static or dynamic.
  bool isStatic;
  /// Pointer to opitional anisotropic mass.
  Vec* anMass;
  /// An optional displacement (useful for parallel quantum wells in 2D).
  double displace;
};

std::ostream& operator<<(std::ostream& os, const Species& s);
#endif
