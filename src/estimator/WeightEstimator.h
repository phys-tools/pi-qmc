#ifndef __WeightEstimator_h_
#define __WeightEstimator_h_
#include "stats/ScalarEstimator.h"
class Paths;
class ScalarAccumulator;

class WeightEstimator: public ScalarEstimator {
public:
    WeightEstimator(ScalarAccumulator*);
    virtual ~WeightEstimator();

    virtual double calcValue();
    virtual void reset();
    virtual void evaluate(const Paths& paths);
};

#endif
