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
#ifndef __RingLattice_h_
#define __RingLattice_h_

#include "Positioner.h"
#include <cstdlib>
#include <blitz/array.h>
class Paths;
class Species;
class MPIManager;

/** Class for positioning particles on a 2-d ring lattice.
  * Require each site is occupied by only 1 particle.
*/
class RingLattice : public Positioner {
public:
  /// Typedefs.
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  RingLattice(Paths& paths, double radius, double angle0, double anglef, 
              double anglex, MPIManager*);
  /// Constructor for placing a single species.
  RingLattice(Paths& paths, double radius, double angle0, double anglef, 
              double anglex, const Species&, MPIManager*);
  ///Virtual destructor
  virtual ~RingLattice() {}
  /// Algorithm run method.
  virtual void run();
private:
  /// Radius of the ring.
  const double radius;
  /// Angle offset
  const double anglex;
  /// The beginning angle
  const double angle0;
  /// The ending angle
  const double anglef;
  /// The index of the first particle.
  const int ifirst;
  /// The number of particles.
  const int npart;
  ///Reference to the Paths.
  Paths& paths;
  /// The MPIManager;
  MPIManager *mpi;
};

#endif
