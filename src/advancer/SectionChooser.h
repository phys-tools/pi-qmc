#ifndef __SectionChooser_h_
#define __SectionChooser_h_

#include "algorithm/CompositeAlgorithm.h"
#include <gsl/gsl_qrng.h>
template<int TDIM> class Beads;
class Paths;
class Permutation;
class Action;
class BeadFactory;

/// Algorithm class for choosing a section.
class SectionChooser: public CompositeAlgorithm {
public:
    SectionChooser(int nlevel, int npart, Paths &paths, Action& action,
            const BeadFactory&);
    virtual ~SectionChooser();

    virtual void run();

    Beads<NDIM>& getBeads() const {
        return *beads;
    }
    Permutation& getPermutation() const {
        return *permutation;
    }
    int getFirstSliceIndex() const {
        return iFirstSlice;
    }
    int getNLevel() const {
        return nlevel;
    }
    const Paths& getPaths() const {
        return *paths;
    }

protected:
    Paths *paths;
    Action *action;
    mutable Beads<NDIM> *beads;
    mutable Permutation *permutation;
    const int nlevel;
    int iFirstSlice;
    gsl_qrng *qrng;
};
#endif
