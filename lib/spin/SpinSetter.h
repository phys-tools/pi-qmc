// $Id: SpinSetter.h,v 1.1 2007/01/03 02:45:49 jshumwa Exp $
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
#ifndef __SpinSetter_h_
#define __SpinSetter_h_

#include "Positioner.h"
#include <blitz/array.h>
class Paths;
class Species;
class MPIManager;

/** Class for initializing spin states of particles.
 * @version $Revision: 1.1 $
 * @author John Shumway */
class SpinSetter : public Positioner {
public:
  /// Typedefs.
  typedef blitz::TinyVector<double,4> SVec;
  /// Constructor.
  SpinSetter(Paths& paths, MPIManager*);
  /// Virtual destructor.
  virtual ~SpinSetter() {}
  /// Algorithm run method.
  virtual void run();
private:
  /// Reference to the Paths.
  Paths& paths;
  /// The number of particles.
  const int npart;
  /// The MPIManager.
  MPIManager *mpi;
};

#endif
