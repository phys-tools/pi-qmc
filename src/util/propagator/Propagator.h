#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

class Propagator {
public:
    Propagator();
    virtual ~Propagator();

    double evaluate();
};

#endif
