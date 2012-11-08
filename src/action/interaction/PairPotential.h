#ifndef __PairPotential_h_
#define __PairPotential_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
class Species;
#include <cstdlib>
#include <blitz/array.h>

/// Class for setting up empirical pair action between particles.
/// We may add methods for parsing formulas later.
/// @author John Shumway. 
class PairPotential {
public:
    virtual ~PairPotential() {}
    virtual double operator()(double r) const=0;
    double getScatteringLength(double mu, double rmax, double dr) const;
};

#endif
