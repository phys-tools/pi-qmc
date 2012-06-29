#ifndef __RandomPermutationChooser_h_
#define __RandomPermutationChooser_h_

#include "advancer/PermutationChooser.h"

/// Class for randomly selecting a permutation.
/// @bug Temporarily set to three particle even permutations.

class RandomPermutationChooser: public PermutationChooser {
public:

    RandomPermutationChooser(const int nsize);

    virtual ~RandomPermutationChooser() {
    }

    virtual bool choosePermutation();
    virtual double getLnTranProb() const {
        return 0.0;
    }
    virtual void init() {
    }
private:
    int nsize;
};
#endif
