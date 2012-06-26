#ifndef __Demo_h_
#define __Demo_h_

#include <iostream>

/// Base class for demonstration input-file generators.
class Demo {
public:
    virtual ~Demo() {}
    /// Construct a demo by name.
    static Demo* getDemo(const std::string& name);
    /// List demos.
    static void listDemos(std::ostream& out);
    /// Method to generate demo.
    virtual void generate() const=0;
};
#endif
