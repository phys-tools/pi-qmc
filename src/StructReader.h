// $Id$
/*  Copyright (C) 2007 John B. Shumway, Jr.

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
#ifndef __StructReader_h_
#define __StructReader_h_

#include "Positioner.h"
#include <blitz/array.h>
class Paths;
class Species;
class MPIManager;

/** Class for position particles from a struct.h5 file.
 * @version $Revision$
 * @author John Shumway */
class StructReader : public Positioner {
public:
  /// Typedefs.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  StructReader(Paths& paths, const std::string &filename,
               MPIManager*);
  /// Virtual destructor.
  virtual ~StructReader() {}
  /// Algorithm run method.
  virtual void run();
private:
  /// Reference to the Paths.
  Paths& paths;
  /// The name of the struct.h5 file.
  const std::string filename;
  /// The MPIManager.
  MPIManager *mpi;
};

#endif
