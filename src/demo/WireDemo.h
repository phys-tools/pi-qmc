#ifndef __WireDemo_h_
#define __WireDemo_h_

#include "Demo.h"
#include <iostream>

class WireDemo: public Demo {
public:
    WireDemo() {
    }
    virtual void generate() const;
};
#endif
