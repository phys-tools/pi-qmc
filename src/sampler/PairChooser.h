#ifndef __PairChooser_h_
#define __PairChooser_h_

#include "ParticleChooser.h"
#include "PermutationChooser.h"
class MultiLevelSampler;
class Species;
class SimulationInfo;
class WalkingChooser;

/// Choose permutations for two species.
class PairChooser: public ParticleChooser, public PermutationChooser {
public:
    PairChooser(const int npart, const Species&, const Species&,
            const int nlevel, const SimulationInfo&);
    virtual ~PairChooser();

    virtual void chooseParticles();

    virtual bool choosePermutation();

    virtual void init();

    virtual double getLnTranProb() const;

    virtual void setMLSampler(const MultiLevelSampler*);
private:
    int npart;
    WalkingChooser &chooser1, &chooser2;
};
#endif
