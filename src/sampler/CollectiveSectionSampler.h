#ifndef __CollectiveSectionSampler_h_
#define __CollectiveSectionSampler_h_
template <int TDIM> class Beads;
class Action;
class BeadFactory;
class AccRejEstimator;
class Paths;
//class SectionChooser;
class SuperCell;
#include "Algorithm.h"
#include "CollectiveSectionMover.h"
#include "SectionSamplerInterface.h"
#include "SectionChooser.h"

class CollectiveSectionSampler : public Algorithm,
                                 public SectionSamplerInterface {
public:
  CollectiveSectionSampler(int, SectionChooser&, Action*, int, 
                      const BeadFactory&, CollectiveSectionMover*, SuperCell*);
  virtual ~CollectiveSectionSampler();

  virtual void run();

  using SectionSamplerInterface::getSectionBeads;
  virtual Beads<NDIM>& getSectionBeads() {return *sectionBeads;}
  virtual const Beads<NDIM>& getSectionBeads() const {return *sectionBeads;}
  virtual Beads<NDIM>& getSectionBeads(int i) {return *sectionBeads;}
  virtual const Beads<NDIM>& getSectionBeads(const int i) const {return *sectionBeads;}

  using SectionSamplerInterface::getMovingBeads;
  virtual Beads<NDIM>& getMovingBeads() {return *movingBeads;}
  virtual const Beads<NDIM>& getMovingBeads() const {return *movingBeads;}
  virtual Beads<NDIM>& getMovingBeads(int i) {return *movingBeads;}
  virtual const Beads<NDIM>& getMovingBeads(const int i) const {return *movingBeads;}

  using SectionSamplerInterface::getMovingIndex;
  virtual const IArray& getMovingIndex() const {return *movingIndex;}
  virtual const IArray& getMovingIndex(const int i) const {return *movingIndex;}

//  virtual int getFirstSliceIndex() const {return 0;}
  virtual const SuperCell& getSuperCell() const {return *cell;}
  virtual bool isSamplingBoth() const {return false;}

  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);

  virtual int getFirstSliceIndex() const {return sectionChooser.getFirstSliceIndex();}


protected:
  bool tryMove();
  /// Number of moving level = log2(nsectionSlice).
  const int nlevel;
  /// Always move all particles.
  const int npart;
  /// Reference to all beads in the section.
  Beads<NDIM> *sectionBeads;
  /// Reference to moving beads.
  Beads<NDIM> *movingBeads;
  /// Index of moving particles = IdentityIndex. Used in action.
  IArray *movingIndex, identityIndex;
  /// Reference to the algorithm that selected the section.
  SectionChooser& sectionChooser;

  SuperCell* cell;
  Action *action;
  const int nrepeat;
  const int maximumMovingCount;
  CollectiveSectionMover* mover;
  AccRejEstimator* accRejEst;
};
#endif
