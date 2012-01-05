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
#ifndef __DoubleSectionChooser_h_
#define __DoubleSectionChooser_h_

#include "SectionChooser.h"
template <int TDIM> class Beads;
class Paths;
class Permutation;
class DoubleAction;
class BeadFactory;

/// Algorithm class for choosing a double section.
/// @version $Revision$
/// @author John Shumway
class DoubleSectionChooser : public SectionChooser {
public:
  /// Constructor.
  DoubleSectionChooser(const int nlevel, Paths&, Action&, DoubleAction&,
    const BeadFactory&);
  /// Virtual destructor.
  virtual ~DoubleSectionChooser();
  /// Run the algorithm.
  virtual void run();
  /// Set the active section for SectionChooser.
  void activateSection(const int i);
  /// Get the section beads.
  Beads<NDIM>& getBeads(const int i) const {return i==1?*beads1:*beads2;}
  /// Get the section permutation.
  Permutation& getPermutation(const int i) const {
    return i==1?*permutation1:*permutation2;}
  /// Get the index of the first slice of the section.
  int getFirstSliceIndex(const int i) const {
    return i==1?iFirstSlice1:iFirstSlice2;}
private:
  /// The action (SectionChooser initializes Action objects for the section).
  DoubleAction& doubleAction;
  /// The section beads.
  mutable Beads<NDIM> *beads1,*beads2; 
  /// The section permutation.
  mutable Permutation *permutation1,*permutation2;
  /// The index of the first slice of the section.
  int iFirstSlice1,iFirstSlice2;
};
#endif
