#ifndef __SHOAction_h_
#define __SHOAction_h_
class SectionSamplerInterface;
class DisplaceMoveSampler;
class Species;
template<int TDIM> class Beads;
#include <cstdlib>
#include <blitz/array.h>
#include "Action.h"

/** Class for calculating the action for a simple harmonic oscillator.
 * The thermal density matrix for a simple harmonic oscillator is,
 * @f[ \rho(x,x';\tau) = \sqrt{\frac{m\omega}{2\pi\sinh\omega\tau}}
 *        \exp\frac{-m\omega[(x^2+{x'}^2)\cosh\omega\tau-2xx']
 *                   }{2\sinh\omega\tau}. @f]
 * For a derivation, see Feynman, <em>Statistical Mecahnics</em> (1972).
 * This class only represents the potential contribution to the
 * action, defined as @f$U=-\log(\rho/\rho^0)@f$.
 * This action is
 * @f[ U(x,x';\tau) = \frac{1}{2}\log\frac{\sinh\omega\tau}{\omega\tau}
 *  +\frac{m\omega[(x^2+{x'}^2)\cosh\omega\tau-2xx']}{2\sinh\omega\tau}
 *  -\frac{m(x-x')^2}{2\tau}. @f]
 * We also need the time derivative,
 * @f[ \dot{U}(x,x';\tau) =
 *   \frac{\omega}{2}\coth\omega\tau
 *  +\frac{m\omega^2}{2}\left[(x^2+{x'}^2)(1-\coth^2\omega\tau)
 *  +\frac{2xx'\cosh\omega\tau}{\sinh^2\omega\tau}\right]
 *  +\frac{m(x-x')^2}{2\tau^2}-\frac{1}{2\tau}. @f]
 *
 * For the virial estimator, we neee to calculate forces.
 * @f[ \nabla_{x}U(x,x';\tau)=
 *   \frac{m\omega[x\cosh\omega\tau-x']}{\sinh\omega\tau}
 *   -\frac{m(x-x')}{\tau} @f]
 * (These formulas generalize to several dimensions by using vector
 * notation.)
 * @todo Treat masses correctly.
 * @todo Calculate forces for virial estimator.
 * @bug Doesn't use mass correctly.
 * @author John Shumway. */
class SHOAction: public Action {
public:
    typedef blitz::Array<int, 1> IArray;
    SHOAction(const double tau, const double omega, const double mass,
            const int ndim, const Species&, const Vec&);
    virtual ~SHOAction() {
    }

    virtual double getActionDifference(const SectionSamplerInterface&,
            int level);
    virtual double getActionDifference(const Paths&,
            const VArray &displacement, int nmoving, const IArray &movingIndex,
            int iFirstSlice, int iLastSlice);
    virtual double getTotalAction(const Paths&, const int level) const;
    virtual void getBeadAction(const Paths&, int ipart, int islice, double& u,
            double& utau, double& ulambda, Vec& fm, Vec& fp) const;
private:
    const double tau;
    const double omega;
    const double mass;
    /// The number of dimensions (can be less than NDIM, i.e., to make a wire).
    const int ndim;
    const int ifirst;
    const int npart;
    /// The location of the minimum of the potential.
    const Vec center;
};
#endif
