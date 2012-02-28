#ifndef __DoubleCollectiveSectionSampler_h_
#define __DoubleCollectiveSectionSampler_h_
template <int TDIM> class Beads;
class DoubleAction;
class Action;
class BeadFactory;
class Paths;
class DoubleSectionChooser;
class CollectiveSectionMover;
class SuperCell;
#include "SectionSamplerInterface.h"
#include "CollectiveSectionSampler.h"

class DoubleCollectiveSectionSampler : public CollectiveSectionSampler {
public:
  DoubleCollectiveSectionSampler(int, DoubleSectionChooser&,
		  Action*, DoubleAction*, int nrepeat, const BeadFactory&,
		  CollectiveSectionMover*, bool, SuperCell*);
  virtual ~DoubleCollectiveSectionSampler();

  virtual void run();

  using CollectiveSectionSampler::getMovingBeads;
  virtual Beads<NDIM>& getMovingBeads(const int i) {
    return i==1?*movingBeads1:*movingBeads2;}
  virtual const Beads<NDIM>& getMovingBeads(const int i) const {
    return i==1?*movingBeads1:*movingBeads2;}

  using CollectiveSectionSampler::getSectionBeads;
  virtual Beads<NDIM>& getSectionBeads(const int i) {
    return i==1?*sectionBeads1:*sectionBeads2;}
  virtual const Beads<NDIM>& getSectionBeads(const int i) const {
    return i==1?*sectionBeads1:*sectionBeads2;}

  using CollectiveSectionSampler::getMovingIndex;
  virtual IArray&  getMovingIndex(const int i) {
    return i==1?*movingIndex1:*movingIndex2;}
  virtual const IArray&  getMovingIndex(const int i) const {
    return i==1?*movingIndex1:*movingIndex2;}
  
  virtual bool isSamplingBoth() const {return samplingBoth;}

  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);

private:
  bool tryMove();
  Beads<NDIM> *sectionBeads1, *sectionBeads2;
  Beads<NDIM> *movingBeads1, *movingBeads2;
  IArray *movingIndex1, *movingIndex2;
  DoubleAction *doubleAction;
  bool samplingBoth;
};
#endif
