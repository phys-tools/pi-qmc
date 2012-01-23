#ifndef __DoubleCollectiveSectionSampler_h_
#define __DoubleCollectiveSectionSampler_h_
template <int TDIM> class Beads;
class DoubleAction;
class Action;
class BeadFactory;
class Paths;
class DoubleSectionChooser;
class CollectiveSectionMover;
#include "SectionSamplerInterface.h"

class DoubleCollectiveSectionSampler : public SectionSamplerInterface {
public:
  DoubleCollectiveSectionSampler(int maximumMovingCount,
          Paths*, DoubleSectionChooser*,
		  Action*, DoubleAction*, const BeadFactory*);
  virtual ~DoubleCollectiveSectionSampler();

  virtual void run();

  virtual const Beads<NDIM>& getSectionBeads() const;
  virtual const Beads<NDIM>& getMovingBeads() const;
  virtual Beads<NDIM>& getSectionBeads();
  virtual Beads<NDIM>& getMovingBeads();
  virtual const IArray& getMovingIndex() const;
  virtual int getFirstSliceIndex() const;
  virtual const SuperCell& getSuperCell() const;


private:
  const int maximumMovingCount;
  CollectiveSectionMover* mover;
};
#endif
