#ifndef __PeriodicGaussian_h_
#define __PeriodicGaussian_h_
#include <cstdlib>
#include <cmath>

/** Function object for a periodic Gaussian.
 * Call evalauate first, the use the getter methods for the function
 * value and derivatives.
 *
 * We start from the k-space representation of a Gaussian,
 * @f[ \exp(-\alpha x^2) =
 * \frac{1}{\sqrt{4\pi\alpha}} \int_{-\infty}^\infty
 * \exp\left(-\frac{k^2}{4 \alpha}\right) \exp(ikx)\, dk. @f]
 * To make this periodic function f(x) on the interval
 * @f$ -L/2 < f < L/2 @f$, we restrict k to integer values,
 * @f$ k_n = 2\pi n/L @f$, where n runs over all integers.
 * Now the integral is a series,
 * @f[ f(x) = \sqrt{\frac{\pi^2}{\alpha L^2}} \left[
 * 1 + \Re e\, 2\sum_{n=1}^\infty
 * \exp\left(-\frac{\pi^2n^2}{\alpha L^2}\right)
 * \exp\left( \frac{2\pi i n x}{L^2} \right) \right ]. @f]
 * The quantity in the square brackets in a theta function
 * with real value `nome'
 * @f$ q = \exp\left(-\frac{\pi^2n^2}{\alpha L^2}\right) @f$
 * and complex argument
 * @f$ w = \exp\left( \frac{2\pi i n x}{L^2} \right) @f$,
 * @f[ f(x) = \sqrt{\frac{\pi^2}{\alpha L^2}} \left[
 * 1 + \Re e\, 2\sum_{n=1}^\infty q^{n^2} w^n \right] @f]
 * The series converges quickly in n, and we choose
 * @f$ n_{max} \ge 2\sqrt{\alpha L^2} @f$, so that
 * @f$ q^{n_{max}^2} < \exp(-4\pi^2) < 10^{-17}@f$.
 * Since we are working in k-space, the first and second derivatives
 * are easily found by multiplying each term by @f$ ik @f$ and
 * @f$ -k^2 @f$, respectively.
 **/
class PeriodicGaussian {
public:
    PeriodicGaussian(double alpha, double length);
    static int numberOfTerms(double alpha, double length);
    double evaluate(double x) const;
    double getValue() const {
        return value;
    }
    double getGradient() const {
        return gradient;
    }
    double getSecondDerivative( ) const {
        return secondDerivative;
    }

private:
    mutable double value;
    mutable double gradient;
    mutable double secondDerivative;
    const double alpha;
    const double length;
    const double nome;
    const double prefactor;
    const double k;
    const int nmax;
    static const double PI;
};
#endif
