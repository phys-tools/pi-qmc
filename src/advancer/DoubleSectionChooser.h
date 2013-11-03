#ifndef __DoubleSectionChooser_h_
#define __DoubleSectionChooser_h_

#include "SectionChooser.h"
template<int TDIM> class Beads;
class Paths;
class Permutation;
class DoubleAction;
class BeadFactory;

/// Algorithm class for choosing a double section.
/// @version $Revision$
/// @author John Shumway
class DoubleSectionChooser: public SectionChooser {
public:
    /// Constructor.
    DoubleSectionChooser(int nlevel, int npart, Paths&, Action&, DoubleAction&,
            const BeadFactory&);
    /// Virtual destructor.
    virtual ~DoubleSectionChooser();
    /// Run the algorithm.
    virtual void run();
    /// Set the active section for SectionChooser.
    void activateSection(const int i);
    /// Get the section beads.
    Beads<NDIM>& getBeads(const int i) const {
        return i == 1 ? *beads1 : *beads2;
    }
    /// Get the section permutation.
    Permutation& getPermutation(const int i) const {
        return i == 1 ? *permutation1 : *permutation2;
    }
    /// Get the index of the first slice of the section.
    int getFirstSliceIndex(const int i) const {
        return i == 1 ? iFirstSlice1 : iFirstSlice2;
    }
private:
    /// The action (SectionChooser initializes Action objects for the section).
    DoubleAction& doubleAction;
    /// The section beads.
    mutable Beads<NDIM> *beads1, *beads2;
    /// The section permutation.
    mutable Permutation *permutation1, *permutation2;
    /// The index of the first slice of the section.
    int iFirstSlice1, iFirstSlice2;
};
#endif
