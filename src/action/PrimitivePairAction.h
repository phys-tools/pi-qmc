#ifndef __EmpriricalInteraction_h_
#define __EmpriricalInteraction_h_
class SectionSamplerInterface;
class Paths;
class Species;
class SimulationInfo;
class PairPotential;
#include "PairAction.h"
#include <cstdlib>
#include <blitz/array.h>

/** Class for setting up empirical pair action between particles.
 * Right now it just uses the primitive approximation.
 * We may add methods for squaring or parsing formulas later. */
class PrimitivePairAction: public PairAction::EmpiricalPairAction {
public:
    virtual double u(double r, int iorder) const;
    virtual double utau(double r, int iorder) const;

    //Construct by providing an empirical potential and time step.
    PrimitivePairAction(const PairPotential& v, const double tau);
    virtual ~PrimitivePairAction() {}

    //Calculate the scattering length.
    double getScatteringLength(Species s1, Species s2, double rmax,
            double dr) const;
private:
    const PairPotential& v;
    const double tau;
};
#endif
