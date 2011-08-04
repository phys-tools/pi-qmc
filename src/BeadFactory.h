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
#ifndef __BeadFactory_h_
#define __BeadFactory_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <iostream>
#include <vector>
#include <string>
#include "Beads.h"

/// Factory for Beads.
/// Uses cloneBase method of Beads to produce specialized Beads for the 
/// simulation.
/// @version $Revision$
/// @author John Shumway 
class BeadFactory {
public:
  /// Constructor.
  BeadFactory();
  /// Copy constructor.
  BeadFactory(const BeadFactory&);
  /// Destructor
  ~BeadFactory();
  /// Get new beads.
  Beads<NDIM>* getNewBeads(const int npart, const int nslice) const;
  /// Add prototype beads to factory.
  int addAuxBeads(BeadsBase*, const std::string& name);
private:
  /// Array of bead names.
  std::vector<std::string> beadNames;
  /// Prototype beads.
  Beads<NDIM>& protoBeads;
};
#endif
