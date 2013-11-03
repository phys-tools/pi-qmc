#ifndef __PermutationChooser_h_
#define __PermutationChooser_h_
class Permutation;
class MultiLevelSampler;

/// Base class for algorithms for selecting a permutation.
/// Defaults to the trivial choice of the direct permutation.
class PermutationChooser {
public:
    PermutationChooser(int size);
    virtual ~PermutationChooser();

    virtual bool choosePermutation();

    const Permutation& getPermutation() {
        return *permutation;
    }

    virtual void init() {
    }

    virtual void setMLSampler(const MultiLevelSampler*) {
    }

    virtual double getLnTranProb() const {
        return 0.0;
    }
protected:
    Permutation* permutation;
};
#endif
