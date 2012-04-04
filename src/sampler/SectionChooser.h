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
#ifndef __SectionChooser_h_
#define __SectionChooser_h_

#include "algorithm/CompositeAlgorithm.h"
#include <gsl/gsl_qrng.h>
template <int TDIM> class Beads;
class Paths;
class Permutation;
class Action;
class BeadFactory;

/// Algorithm class for choosing a section.
/// @version $Revision$
/// @author John Shumway
class SectionChooser : public CompositeAlgorithm {
public:
  /// Constructor.
  SectionChooser(int nlevel, int npart, Paths &paths, Action& action,
                 const BeadFactory&);
  /// Virtual destructor.
  virtual ~SectionChooser();
  /// Run the algorithm.
  virtual void run();
  /// Get the section beads.
  Beads<NDIM>& getBeads() const {return *beads;}
  /// Get the section permutation.
  Permutation& getPermutation() const {return *permutation;}
  /// Get the index of the first slice of the section.
  int getFirstSliceIndex() const {return iFirstSlice;}
  /// Get the number of levels in the section.
  int getNLevel() const {return nlevel;}
  /// Get a reference to the paths.
  const Paths& getPaths() const {return *paths;}
protected:
  /// The paths.
  Paths *paths;
  /// The action (SectionChooser initializes Action objects for the section).
  Action *action;
  /// The section beads.
  mutable Beads<NDIM> *beads;
  /// The section permutation.
  mutable Permutation *permutation;
  /// The number of levels.
  const int nlevel;
  /// The index of the first slice of the section.
  int iFirstSlice;
  /// The quasirandom number generator.
  gsl_qrng *qrng;
};
#endif
