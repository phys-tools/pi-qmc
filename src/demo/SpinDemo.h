#ifndef __SpinDemo_h_
#define __SpinDemo_h_

#include "Demo.h"
#include <iostream>

/// Class for single spin demo.
class SpinDemo: public Demo {
public:
    virtual ~SpinDemo() {}
    virtual void generate() const;
};
#endif
