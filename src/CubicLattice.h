// $Id: CubicLattice.h,v 1.7 2007/11/26 00:00:14 jshumwa Exp $
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
#ifndef __CubicLattice_h_
#define __CubicLattice_h_

#include "Positioner.h"
#include <blitz/array.h>
class Paths;
class Species;
class MPIManager;

/** Class for positioning particles in a cubic lattice.
 * @version $Revision: 1.7 $
 * @author John Shumway */
class CubicLattice : public Positioner {
public:
  /// Typedefs.
  typedef blitz::TinyVector<int,NDIM> IVec;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  CubicLattice(Paths& paths, double a, double scatter,
               const IVec nmax, const Vec r0, MPIManager*);
  /// Constructor for placing a single species.
  CubicLattice(Paths& paths, double a, double scatter,
               const IVec nmax, const Vec r0, 
               const Species&, MPIManager*);
  /// Virtual destructor.
  virtual ~CubicLattice() {}
  /// Algorithm run method.
  virtual void run();
private:
  /// Reference to the Paths.
  Paths& paths;
  /// Lattice constant.
  const double a;
  /// Percentage random displacement.
  const double scatter;
  /// Number of grid points in each direction.
  const IVec nmax;
  /// Center of grid.
  const Vec r0;
  /// The index of the first particle.
  const int ifirst;
  /// The number of particles.
  const int npart;
  /// The MPIManager.
  MPIManager *mpi;
};

#endif
