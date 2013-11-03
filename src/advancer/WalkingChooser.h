#ifndef __WalkingChooser_h_
#define __WalkingChooser_h_

#include "PermutationChooser.h"
#include "SpeciesParticleChooser.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
class MultiLevelSampler;
class SimulationInfo;
class PeriodicGaussian;

class WalkingChooser:
    public PermutationChooser, public SpeciesParticleChooser {
public:
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::Array<double, 2> Mat;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<PeriodicGaussian*, 1> PGArray;
    WalkingChooser(const int nsize, const Species&, const int nlevel,
            const SimulationInfo&);
    virtual ~WalkingChooser();
    virtual void setMLSampler(const MultiLevelSampler*);
    virtual bool choosePermutation();
    virtual void chooseParticles();
    virtual void init();
    int iSearch(int part, double x);
    virtual double getLnTranProb() const {
        return log(prob);
    }
protected:
    Mat t;
    Mat cump;
    const MultiLevelSampler *multiLevelSampler;
    PGArray pg;
private:
    int nsize;
    const double mass;
    const double tau;
    double prob;
};
#endif
