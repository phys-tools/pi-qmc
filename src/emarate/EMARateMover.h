#ifndef __EMARateMover_h_
#define __EMARateMover_h_
#include <cstdlib>
#include <blitz/array.h>
#include <vector>
#include "sampler/Mover.h"
#include "sampler/PermutationChooser.h"
#include "sampler/ParticleChooser.h"
class SimulationInfo;
class PeriodicGaussian;

class EMARateMover : public Mover,
                     public PermutationChooser,
                     public ParticleChooser {
public:
    typedef blitz::Array<int,1> IArray;
    typedef blitz::Array<double,1> Array;
    typedef blitz::Array<PeriodicGaussian*,3> PGArray;
    typedef blitz::TinyVector<double,NDIM> Vec;
    EMARateMover(const SimulationInfo&, const int maxlevel, const double pgDelta);
    virtual ~EMARateMover();
    /// Move the samplers moving beads for a given level, returning
    /// the probability for the old move divided by the probability for the
    /// new move.
    virtual double makeMove(MultiLevelSampler&, const int level);
    virtual double makeDelayedMove(MultiLevelSampler&, const int level) ;
    virtual double getForwardProb() {return forwardProb;}
    virtual void chooseParticles() {}
private:
    /// The inverse mass, @f$\lambda=1/2m@f$.
    blitz::Array<double,1> lambda;
    double tau;
    double forwardProb;
    bool isSamplingRadiating;
    double earlierTransitions;
    const double C;
};
#endif
