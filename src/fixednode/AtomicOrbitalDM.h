#ifndef ATOMICORBITALDM_H_
#define ATOMICORBITALDM_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>

/// Classes for atomic orbital density matricies.
/// These are optimized to work for vectors of identical orbitals.
class AtomicOrbitalDM {
public:
    typedef blitz::TinyVector<double,NDIM> Vec;
    typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;
    typedef blitz::Array<double,2> Matrix;
    typedef blitz::Array<Vec,2> VMatrix;
    typedef blitz::Array<double,1> Array;
    typedef blitz::Array<int,1> IArray;
    typedef blitz::Array<int,2> IArray2;
    typedef blitz::Array<double,2> Array2;
    typedef blitz::Array<Vec,2> VArray2;
    typedef blitz::ColumnMajorArray<2> ColMajor;

  AtomicOrbitalDM(int ifirst, int npart, int nfermion, double weight);
  const double weight;
  const int ifirst;
  const int npart;
  /// Vector displacement of particle i from orbital center j in section 1.
  mutable VArray2 dist1;
  /// Vector displacement of particle i from orbital center j in section 2.
  mutable VArray2 dist2;
  /// Get the array of vector distances for section 1.
  VArray2& getD1Array() const {return dist1;}
  /// Get the array of vector distances for section 2.
  VArray2& getD2Array() const {return dist2;}
  /// Structure for value and gradient.
  struct ValAndGrad {
    ValAndGrad() : val(0.), grad1(0.), grad2(0.) {}
    double val; Vec grad1, grad2;
    ValAndGrad operator+=(ValAndGrad& s) {
      val += s.val; grad1 += s.grad1; grad2 += s.grad2; return *this;
    }
  };
  /// Calculate matrix values using vector distance arrays.
  /// On exit, mat(i,j) has matrix elements for particles ri and rj'.
  virtual void evaluateValue(Matrix&, double scale) const {;}
  /// Calculate matrix values using vector distance arrays.
  /// On exit work and distance arrays are setup for ValAndGrad operator().
  virtual void evaluateValueAndGrad() const {;}
  /// Return value and gradients of matrix element between ri  and rj'.
  /// Must call evaluateValueAndGrad() first to initialize work arrays.
  virtual ValAndGrad operator()(int i, int j) const {
    ValAndGrad temp; return temp;}
};
#endif
