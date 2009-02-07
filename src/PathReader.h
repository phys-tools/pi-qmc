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
#ifndef __PathReader_h_
#define __PathReader_h_

#include <string>
#include "Positioner.h"
class Paths;
class MPIManager;
class BeadFactory;

/** Class for positioning particles in a cubic lattice.
 * @version $Revision$
 * @author John Shumway */
class PathReader : public Positioner {
public:
  /// Constructor.
  PathReader(Paths& paths, std::string& filename, 
             const BeadFactory &beadFactory,
             const int bfactor=1, MPIManager *mpi=0)
   : paths(paths), filename(filename), beadFactory(beadFactory),
     bfactor(bfactor), mpi(mpi) {}
  /// Virtual destructor.
  virtual ~PathReader() {}
  /// Algorithm run method.
  virtual void run();
private:
  /// Reference to the Paths.
  Paths& paths;
  /// Input file name.
  const std::string filename;
  /// Bead factory.
  const BeadFactory& beadFactory;
  /// Optional beta multiplication factor for restarting at lower temperature.
  const int bfactor;
  /// Pointer to the MPI manager. 
  MPIManager *mpi;
};

#endif
