#ifndef ESTIMATORITERATOR_H_
#define ESTIMATORITERATOR_H_

#include <list>
class Estimator;

class EstimatorIterator {
public:
    EstimatorIterator(std::list<Estimator*> &list);
    virtual ~EstimatorIterator();
    bool step();
    Estimator* operator*();
private:
    std::list<Estimator*>::iterator iterator;
    std::list<Estimator*>::iterator end;
};

#endif
