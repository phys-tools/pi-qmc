#ifndef __DoubleCollectiveSectionSampler_h_
#define __DoubleCollectiveSectionSampler_h_
template <int TDIM> class Beads;
class DoubleAction;
class Action;
class BeadFactory;
class Paths;
class DoubleSectionChooser;

class DoubleCollectiveSectionSampler {
public:
  DoubleCollectiveSectionSampler(Paths*, DoubleSectionChooser*,
		  Action*, DoubleAction*, const BeadFactory*);
  virtual ~DoubleCollectiveSectionSampler();

  virtual void run();
private:
};
#endif
