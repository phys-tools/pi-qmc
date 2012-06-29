#ifndef __MiddleSectionChooser_h_
#define __MiddleSectionChooser_h_

#include "SectionChooser.h"

/// Algorithm class for choosing a section.
class MiddleSectionChooser: public SectionChooser {
public:
    MiddleSectionChooser(const int nlevel, Paths &paths, Action& action,
            const BeadFactory&);
    virtual ~MiddleSectionChooser();

    virtual void run();
};
#endif
