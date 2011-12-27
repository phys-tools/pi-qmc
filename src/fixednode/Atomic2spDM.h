#ifndef ATOMIC2SPDM_H_
#define ATOMIC2SPDM_H_

#include "AtomicOrbitalDM.h"

/// Unnormalized density matrix for 2s and 2p orbital.
/// To avoid negative regions, this does not include the radial
/// node in the 2s orbital, and there is enough 2s contribution
/// to avoid an angular node (provided @f$0<p_{\text{angle}}@f<1$).
/// @f[\rho(r,r') = C^2 rr' e^{-Z(r+r')/2}
///                  (1+\hat{\mathbf{r}}\cdot\hat{\mathbf{r}}'),@f]
/// where @f$C=\sqrt{Z^5/32\pi}@f$.
/// and @f$0<p_{\text{angle}}@f<1$ controls the strength
/// angular dependence (@f$p_{\text{angle}}=0@f$ for s-only).
class Atomic2spDM : public AtomicOrbitalDM {
public:
  Atomic2spDM(double Z, int ifirst, int npart, int nfermion,
              double pweight, double weight);
  const int nfermion;
  /// Coefficient before the exponential, @f$ C = \sqrt{Z^5/32\pi} @f$.
  const double coef;
  /// Storage for @f$ \psi_k(r_i) @f$.
  mutable Array2 work1;
  /// Storage for @f$ \psi^*_k(r_i) @f$.
  mutable Array2 work2;
  /// Storage for unit vector r_ik.
  mutable VArray2 work3;
  /// Storage for unit vector r'_ik.
  mutable VArray2 work4;
  virtual void evaluateValue(Matrix&, double scale) const;
  virtual void evaluateValueAndGrad() const;
  virtual ValAndGrad operator()(int i, int j) const;
  const double pweight;
  const double Z;
  static const double PI;
};
#endif
