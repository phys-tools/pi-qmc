#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PairPotential.h"

double PairPotential::getScatteringLength(double mu,
    double rmax, double dr) const {
  //Use Numerov method to integrate SE. a=r-psi(r)/psi'(r).
  double psi = dr, psim=0.;
  double f = 2*mu*operator()(dr), fm = 2*mu*operator()(0.);
  double r = dr;
  while (r<rmax) {
    r += dr;
    double fp = 2*mu*operator()(r);
    double psip = (2.*psi - psim + dr*dr*(fm*psim + 10*f*psi)/12.)
                 /(1.-dr*dr*fp/12.);
    psim=psi; fm=f;
    psi=psip; f=fp;
  }
  return (r-0.5*dr) - 0.5*dr*(psi+psim)/(psi-psim);
}
