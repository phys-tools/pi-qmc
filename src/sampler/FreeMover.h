#ifndef __FreeMover_h_
#define __FreeMover_h_
#include <cstdlib>
#include <blitz/array.h>
#include <vector>
#include "Mover.h"
class SimulationInfo;
class PeriodicGaussian;
/** Select a trial move for beads by free particle sampling.
 * @todo May want to make this a subclass of Action, with a flag
 *       for whether to return probability or allow to cancel
 *       for efficiency.
 * @version $Revision$
 * @author John Shumway. */
class FreeMover: public Mover {
public:
    /// Typedefs.
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::Array<PeriodicGaussian*, 3> PGArray;
    typedef blitz::TinyVector<double, NDIM> Vec;
    /// Construct by providing lambda (@f$\lambda=1/2m@f$) timestep @f$\tau@f$.
    FreeMover(const double lambda, const int npart, const double tau);
    /// Construct by providing simInfo.
    FreeMover(const SimulationInfo&, const int maxlevel, const double pgDelta);
    /// Virtual destructor.
    virtual ~FreeMover();
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
    blitz::Array<double, 1> lambda;
    /// The timestep.
    double tau;
    /// Periodic gaussians.
    PGArray pg;
    IArray specIndex;
    /// forward transition prob
    double forwardProb;
};
#endif
