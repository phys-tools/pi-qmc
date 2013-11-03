#ifndef ATOMIC1SDM_H_
#define ATOMIC1SDM_H_

#include "AtomicOrbitalDM.h"

/// Normalized density matrix for 1s orbital.
/// @f[\rho(r,r') = C^2 e^{-Z(r+r')},@f] where @f$C=\sqrt{Z^3/\pi}@f$.
class Atomic1sDM : public AtomicOrbitalDM {
public:
    Atomic1sDM(double Z, int ifirst, int npart, int nfermion, double weight);
    virtual ~Atomic1sDM();
    const int nfermion;
    /// Coefficient before the exponential, @f$ C = \sqrt{Z^3/\pi} @f$.
    const double coef;
    /// Storage for @f$ \psi_k(r_i) @f$.
    mutable Array2 work1;
    /// Storage for @f$ \psi^*_k(r_i) @f$.
    mutable Array2 work2;
    virtual void evaluateValue(Matrix&, double scale) const;
    virtual void evaluateValueAndGrad() const;
    virtual ValAndGrad operator()(int i, int j) const;
    const double Z;
    static const double PI;
};
#endif
