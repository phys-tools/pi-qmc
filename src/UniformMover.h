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
#ifndef __UniformMover_h_
#define __UniformMover_h_
#include <blitz/array.h>
#include <vector>

class MPIManager;
class DisplaceMoveSampler;
class DoubleDisplaceMoveSampler;

class UniformMover {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;

  UniformMover(double dist, const MPIManager *mpi);
  virtual ~UniformMover();
  virtual double makeMove(VArray&, const int&) const;
  virtual double makeMove(DoubleDisplaceMoveSampler&) const;
private:
  const double dist;
  const MPIManager* mpi;
};
#endif
