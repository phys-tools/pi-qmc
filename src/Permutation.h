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
#ifndef __Permutation_h_
#define __Permutation_h_

#include <blitz/array.h>

/// Class for representing a permuation.
/// @version $Revision$
/// @author John Shumway
class Permutation {
public:
  /// Constructor, defaults to the identity.
  Permutation(int npart);
  /// Copy constructor.
  Permutation(const Permutation&);
  /// Assignment operator.
  Permutation& operator=(const Permutation&);
  /// Return a permuted index.
  int operator[](const int i) const {return permutation(i);}
  /// Return a permuted index.
  int& operator[](const int i) {return permutation(i);}
  /// Reset the permutation to the identity.
  void reset();
  /// Prepend a Permutation to this Permutation.
  Permutation& prepend(const Permutation&);
  /// Apend a Permutation to this Permutation.
  Permutation& append(const Permutation&);
  /// Set to the inverse of another permutation.
  void setToInverse(const Permutation&);
  /// Check if a permutatino is the identity.
  bool isIdentity() const;
  friend std::ostream& operator<<(std::ostream&, const Permutation&);
private:
  /// The permutation.
  blitz::Array<int,1> permutation;
};
#endif
