#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PrimitivePairAction.h"
#include "PairPotential.h"
#include "base/Species.h"

double PrimitivePairAction::u(double r, int iorder) const {
    return (iorder == 0) ? v(r) * tau : 0.;
}

double PrimitivePairAction::utau(double r, int iorder) const {
    return (iorder == 0) ? v(r) : 0.;
}

PrimitivePairAction::PrimitivePairAction(const PairPotential& v,
        const double tau) :
        v(v), tau(tau) {
}

double PrimitivePairAction::getScatteringLength(Species s1, Species s2,
        double rmax, double dr) const {
    //Use Numerov method to integrate SE. a=r-psi(r)/psi'(r).
    double mu = 1. / (1. / s1.mass + 1. / s2.mass);
    double psi = dr, psim = 0.;
    double f = 2 * mu * v(dr), fm = 2 * mu * v(0.);
    double r = dr;
    while (r < rmax) {
        r += dr;
        double fp = 2 * mu * v(r);
        double psip = (2. * psi - psim
                + dr * dr * (fm * psim + 10 * f * psi) / 12.)
                / (1. - dr * dr * fp / 12.);
        psim = psi;
        fm = f;
        psi = psip;
        f = fp;
    }
    return (r - 0.5 * dr) - 0.5 * dr * (psi + psim) / (psi - psim);
}
