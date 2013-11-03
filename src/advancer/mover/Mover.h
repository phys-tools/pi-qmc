#ifndef __Mover_h_
#define __Mover_h_

class MultiLevelSampler;
/// Virtual base class for routines to select trial moves for beads.

class Mover {
public:
    virtual ~Mover() {
    }
    /// Move the samplers moving beads for a given level, returning
    /// the probability for the old move divided by the probability for the
    /// new move.
    virtual double makeMove(MultiLevelSampler&, const int level)=0;
    virtual double makeDelayedMove(MultiLevelSampler&, const int level)=0;
    virtual double getForwardProb() = 0;
};
#endif
