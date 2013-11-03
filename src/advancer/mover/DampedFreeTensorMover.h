#ifndef __DampedFreeTensorMover_h_
#define __DampedFreeTensorMover_h_
#include "Mover.h"
class SimulationInfo;
#include <cstdlib>
#include <blitz/array.h>
/** Select a trial move for beads by free particle sampling.
 * @todo May want to make this a subclass of Action, with a flag
 *       for whether to return probability or allow to cancel
 *       for efficiency.
 * @bug  Needs to calculate the transiton probability correctly when
 *       periodic images are significant.  Currently the routine
 *       neglects this.
 * @bug  Hard-coded for Si/Ge electron and hole.
 * @version $Revision$
 * @author John Shumway. */
class DampedFreeTensorMover: public Mover {
public:
    typedef blitz::TinyVector<double, NDIM> Vec;
    /// Construct by providing lambda (@f$\lambda=1/2m@f$) timestep @f$\tau@f$.
    DampedFreeTensorMover(const SimulationInfo &simInfo,
            const int saturationLevel);
    /// Virtual destructor.
    virtual ~DampedFreeTensorMover() {
    }
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
    /// The inverse mass, @f$\lambda=1/2m@f$.
    blitz::Array<blitz::TinyVector<double, NDIM>, 1> lambda;
    /// The timestep.
    double tau;
    /// The level at which lambda saturation occurs.
    int saturationLevel;
};
#endif
