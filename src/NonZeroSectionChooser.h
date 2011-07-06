// $Id: NonZeroSectionChooser.h 338 2010-11-30 18:56:16Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr.

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
#ifndef __NonZeroSectionChooser_h_
#define __NonZeroSectionChooser_h_

#include "SectionChooser.h"
template <int TDIM> class Beads;
class Paths;
class Permutation;
class Action;
class BeadFactory;

/// Algorithm class for choosing a section.
/// @version $Revision: 338 $
/// @author John Shumway
class NonZeroSectionChooser : public SectionChooser {
public:
  /// Constructor.
  NonZeroSectionChooser(const int nlevel, Paths &paths, Action& action,
                 const BeadFactory&);
  /// Virtual destructor.
  virtual ~NonZeroSectionChooser();
  /// Run the algorithm.
  virtual void run();
};
#endif
