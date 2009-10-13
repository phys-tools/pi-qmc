// $Id$
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#ifndef __DoubleDisplaceMoveSampler_h_
#define __DoubleDisplaceMoveSampler_h_
class DoubleAction;

#include "DisplaceMoveSampler.h"
#include <vector>
#include <blitz/array.h>
#include <iostream>

/** Class to perform classical displacements of the particles.
  @author Saad Khairallah, John Shumway
*/
class DoubleDisplaceMoveSampler : public DisplaceMoveSampler {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  DoubleDisplaceMoveSampler(const int nmoving, const int nrepeat,
    Paths&, ParticleChooser&, const UniformMover&, Action*, DoubleAction*,
    const MPIManager* mpi);
  /// Destructor.
  virtual ~DoubleDisplaceMoveSampler();
protected:
  /// First slice of second beads.
  int iFirstSlice2;
  /// Pointer to the double action.
  DoubleAction* doubleAction;
  /// Method to atempt a Monte Carlo move, return true if accepted.
  virtual bool tryMove();
};
#endif
