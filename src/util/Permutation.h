#ifndef __Permutation_h_
#define __Permutation_h_

#include <cstdlib>
#include <blitz/array.h>

/// Class for representing a permuation.
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
  /// Check if a permutation is the identity.
  bool isIdentity() const;
  /// Calculate the inversion count of the permutation.
  int inversionCount() const;
  /// The sign of the permutation.
  int sign() const;
  friend std::ostream& operator<<(std::ostream&, const Permutation&);
private:
  /// The permutation.
  blitz::Array<int,1> permutation;
};
#endif
