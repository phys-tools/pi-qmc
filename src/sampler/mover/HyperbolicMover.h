#ifndef __HyperbolicMover_h_
#define __HyperbolicMover_h_
#include <cstdlib>
#include <blitz/array.h>
#include <vector>
#include "sampler/mover/Mover.h"
class SimulationInfo;
class PeriodicGaussian;
/** Select a trial move for beads by free particle sampling.
 * @todo May want to make this a subclass of Action, with a flag
 *       for whether to return probability or allow to cancel
 *       for efficiency*/
class HyperbolicMover: public Mover {
public:
    /// Typedefs.
    typedef blitz::Array<double, 1> Array;
    typedef blitz::Array<double, 2> Array2;
    typedef blitz::TinyVector<double, NDIM> Vec;
    /// Construct by providing simInfo.
    HyperbolicMover(const SimulationInfo&, const int maxlevel);
    /// Virtual destructor.
    virtual ~HyperbolicMover();
    /// Move the samplers moving beads for a given level, returning
    /// the probability for the old move divided by the probability for the
    /// new move.
    virtual double makeMove(MultiLevelSampler&, const int level);
    virtual double makeDelayedMove(MultiLevelSampler&, const int level) {
        return 0;
    }
    virtual double getForwardProb() {
        return 0;
    }
private:
    /// The timestep.
    const double tau;
    /// The mass.
    const double mass;
    /// The anisotropy parameter.
    const double alpha;
    /// Sampling scale function lookup table.
    Array2 fGtable, fG2table;
    /// Inverse of grid spacing for scale function table.
    double drinv;
    /// Sampling scale function for g.
    double fG(const double r, const int ilevel) const;
    /// Sampling scale function for g2.
    double fG2(const double r, const int ilevel) const;
    /// Propagator.
    double g(const double r, const double t, const int ilevel) const;
    /// Unscaled gaussian propogator.
    double g0(const double r, const double t) const;
    /// Set scale function table for sampling.
    void setFTable(int maxlevel);
    /// Pi.
    static const double pi;
    /// Norm for g and g2.
    Array normG, normG2;
    /// Sample the displacement from a midpoint.
    Vec sampleDelta(Vec& deltaOver2, const double teff, const int ilevel);
    /// Calculate the sampling probability.
    double g1g2(const Vec center, const Vec delta, const double teff,
            const int level);
};
#endif
