#ifndef __SHODemo_h_
#define __SHODemo_h_

#include "Demo.h"
#include <iostream>

/// Class for simple harmonic oscillator demo.
class SHODemo: public Demo {
public:
    virtual ~SHODemo() {}
    virtual void generate() const;
};
#endif
