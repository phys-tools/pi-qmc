#ifndef __NonZeroSectionChooser_h_
#define __NonZeroSectionChooser_h_

#include "SectionChooser.h"
template<int TDIM> class Beads;
class Paths;
class Permutation;
class Action;
class BeadFactory;

/// Algorithm class for choosing a section.
class NonZeroSectionChooser: public SectionChooser {
public:
    NonZeroSectionChooser(const int nlevel, Paths &paths, Action& action,
            const BeadFactory&);

    virtual ~NonZeroSectionChooser();

    virtual void run();
};
#endif
