#ifndef __SuperCell_h_
#define __SuperCell_h_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <blitz/tinyvec.h>
#include <iostream>

/// A simple rectangular supercell.
/// @version $Revision$
/// @author John Shumway
class SuperCell {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Constructor 
  SuperCell(const Vec a);
  /// Virtual destructor.
  virtual ~SuperCell();
  /// Compute the reciprical lattice vectors.
  void computeRecipricalVectors();
  /// Set to the smallest displacement with PBC.
  /// @todo May want to improve algorithm here.
  Vec& pbc(Vec&) const;
  /// smallest displacement along an axis with PBC
  double pbc(double dist, int idim) const;
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
