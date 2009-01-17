// $Id: SuperCell.h,v 1.4 2006/10/18 17:08:19 jshumwa Exp $
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
#ifndef __SuperCell_h_
#define __SuperCell_h_

#include <blitz/tinyvec.h>
#include <iostream>

/// A simple rectangular supercell.
/// @version $Revision: 1.4 $
/// @author John Shumway
class SuperCell {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor 
  SuperCell(const Vec a);
  /// Virtual destructor.
  virtual ~SuperCell() {}
  /// Compute the reciprical lattice vectors.
  void computeRecipricalVectors();
  /// Set to the smallest displacement with PBC.
  /// @todo May want to improve algorithm here.
  Vec& pbc(Vec&) const;
  /// Write info to an ostream.
  virtual std::ostream& write(std::ostream&) const;
  /// Supercell dimensions.
  const Vec a;
  /// Inverse dimensions.
  Vec b;
  /// Cutoff for projecting to smallest vector in pbc method.
  double rcut2;
  /// Indexed access to box lengths.
  double operator[](const int i) const {return a[i];} 
};

std::ostream& operator<<(std::ostream&, const SuperCell&);
#endif
