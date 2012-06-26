#ifndef __HAtomDemo_h_
#define __HAtomDemo_h_

#include "Demo.h"
#include <iostream>

class HAtomDemo : public Demo {
public:
    virtual ~HAtomDemo() {}
    virtual void generate() const;
};
#endif
