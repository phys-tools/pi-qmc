#ifndef __CaoBerneAction_h_
#define __CaoBerneAction_h_
class SectionSamplerInterface;
class Paths;
class Species;
class SimulationInfo;
class PairPotential;
#include "PairAction.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for setting up Cao-Berne approximation to the hard-sphere action.
 * @f[
 * u(r,r',\cos \theta;\tau)
 * = -\ln\left[
 * 1 - \frac{a(r+r')-a^2}{rr'}
 * e^{-\frac{\mu}{\hbar^2\tau}[rr'+a^2-a(r+r')](1+\cos\theta)}
 *  \right]
 * @f] */
class CaoBerneAction: public PairAction::EmpiricalPairAction {
public:
    typedef blitz::Array<double, 1> Array;
    typedef blitz::Array<double, 2> Array2;
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::TinyVector<double, NDIM> Vec;

    //Construct by providing a reduced mass, radius, time step, and order.
    CaoBerneAction(double mu, double tau, double radius, int norder);
    virtual ~CaoBerneAction() {}

    //Calculate the action.
    double getAction(double r1, double r2, double costheta, double dtau) const;
    void calcU(double q) const;
    virtual double u(double r, int idata) const;
    virtual double utau(double r, int idata) const;
private:
    double radius;
    const double mu;
    const double tau;
    const int norder;
    const int ndata;
    mutable double q;
    mutable Array2 uarray;
    mutable Array2 mat;
    mutable IArray ipiv;
    mutable Array s2, z;
};
#endif
