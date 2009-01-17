// $Id: PermutationChooser.h,v 1.8 2008/11/25 13:11:21 jshumwa Exp $
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
#ifndef __PermutationChooser_h_
#define __PermutationChooser_h_
class Permutation;
/// Base class for algorithms for selecting a permutation.
/// Defaults to the trivial choice of the direct permutation.
/// @version $Revision: 1.8 $
/// @author John Shumway
class PermutationChooser {
public:
  /// Construct by giving the size of the permutation.
  PermutationChooser(int size);
  /// Virtual destructor.
  virtual ~PermutationChooser();
  /// Choose and return acceptance outcome.
  virtual bool choosePermutation()=0;
  /// Return the choosen Permutation.
  const Permutation& getPermutation() {return *permutation;}
  /// Initialize.
  virtual void init() =0;
  virtual double getLnTranProb() const =0;
protected:
  /// The current permutation.
  Permutation* permutation;
};
#endif
