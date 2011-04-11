// $Id: CollectiveMover.h $
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
#ifndef __CollectiveMover_h_
#define __CollectiveMover_h_
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <vector>
#include "UniformMover.h"

class Paths;

class CollectiveMover : public UniformMover {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;

  CollectiveMover(Paths& paths, const int& iFirstSlice, const Vec kvec,
		  const Vec amp, const MPIManager *mpi);
  virtual ~CollectiveMover();
  virtual double makeMove(VArray&, const IArray&) const;
private:
  Paths& paths;
  const int iFirstSlice;
  const Vec kvec, amp;
};
#endif
