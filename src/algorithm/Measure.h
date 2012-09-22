#ifndef __Measure_h_
#define __Measure_h_

#include "Algorithm.h"
class Paths;
class Estimator;
class PartitionWeight;
#include <vector>

/** Algorithm class for measure estimators.
 * @author John Shumway */
class Measure : public Algorithm {
public:
  Measure(Paths&, std::vector<Estimator*>, PartitionWeight* weight);
  virtual ~Measure() {}

  virtual void run();
private:
  Paths& paths;
  const std::vector<Estimator*> estimator;
  PartitionWeight* weight;
};
#endif
