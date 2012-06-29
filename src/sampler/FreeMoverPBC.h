#ifndef __FreeMoverPBC_h_
#define __FreeMoverPBC_h_
#include <cstdlib>
#include <blitz/array.h>
#include <vector>
#include "Mover.h"
class SimulationInfo;
class PeriodicGaussian;
/** Select a trial move for beads by free particle sampling.
 * This class treats periodic boundary conditions more carefully, which
 * is necessary for very small box sizes.
 *
 * The difficulty, is that the bisection algorithm has 2**NDIM midpoints
 * when you consider periodic boundary conditions. If these are not sampled
 * carefully, some windings or permutations may be missed, raising the energy.
 *
 * @version $Revision$
 * @bug Delayed rejection move not properly coded.
 * @author John Shumway. */
class FreeMoverPBC: public Mover {
public:
    /// Typedefs.
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::Array<PeriodicGaussian*, 3> PGArray;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::TinyVector<int, NDIM> IVec;
    /// Construct by providing lambda (@f$\lambda=1/2m@f$) timestep @f$\tau@f$.
    FreeMoverPBC(const double lambda, const int npart, const double tau);
    /// Construct by providing simInfo.
    FreeMoverPBC(const SimulationInfo&, const int maxlevel,
            const double pgDelta);
    /// Virtual destructor.
    virtual ~FreeMoverPBC();
    /// Move the samplers moving beads for a given level, returning
    /// the probability for the old move divided by the probability for the
    /// new move.
    virtual double makeMove(MultiLevelSampler&, const int level);
    virtual double makeDelayedMove(MultiLevelSampler&, const int level);
    virtual double getForwardProb() {
        return forwardProb;
    }
private:
    /// The inverse mass, @f$\lambda=1/2m@f$.
    Array lambda;
    /// The timestep.
    double tau;
    /// Periodic gaussians.
    PGArray pg;
    /// Vector with supercell dimensions.
    Vec length;
    /// Index of species types.
    IArray specIndex;
    /// forward transition prob
    double forwardProb;
};
#endif
