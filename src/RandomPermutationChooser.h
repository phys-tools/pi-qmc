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
#ifndef __RandomPermutationChooser_h_
#define __RandomPermutationChooser_h_

#include "PermutationChooser.h"
/* Class for randomly selecting a permutation.
 * @bug Temporarily set to three particle even permutations.
 * @version $Revision$
 * @author John Shumway */
class RandomPermutationChooser : public PermutationChooser {
public:
  /// Construct by giving the size.
  RandomPermutationChooser(const int nsize);
  /// Virtual destructor.
  virtual ~RandomPermutationChooser() {}
  /// Choose and return acceptance outcome.
  virtual bool choosePermutation();
  virtual double getLnTranProb() const {return 0.0;}
  virtual void init() {}
private:
  /// Size of the permutation.
  int nsize;
};
#endif
