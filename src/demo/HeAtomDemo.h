#ifndef __HeAtomDemo_h_
#define __HeAtomDemo_h_

#include "Demo.h"
#include <iostream>

class HeAtomDemo : public Demo {
public:
    virtual ~HeAtomDemo() {}
    virtual void generate() const;
};
#endif
