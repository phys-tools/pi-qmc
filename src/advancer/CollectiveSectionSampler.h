#ifndef __CollectiveSectionSampler_h_
#define __CollectiveSectionSampler_h_
template<int TDIM> class Beads;
class Action;
class BeadFactory;
class AccRejEstimator;
class Paths;
class SuperCell;
#include "algorithm/Algorithm.h"
#include "CollectiveSectionMover.h"
#include "SectionSamplerInterface.h"
#include "SectionChooser.h"

/** Sampling a segment of paths with a specified mover, aiming to improve
 the efficiency for dynamic estimators. */
class CollectiveSectionSampler: public Algorithm, public SectionSamplerInterface {
public:
    CollectiveSectionSampler(int, SectionChooser&, Action*, int,
            const BeadFactory&, CollectiveSectionMover*, SuperCell*);
    virtual ~CollectiveSectionSampler();

    virtual void run();
    bool tryMove();
    void initializeMovingBeads();
    void reportAtempt() const;
    void reportAcceptance() const;
    void updateSectionBeads();

    using SectionSamplerInterface::getSectionBeads;
    virtual Beads<NDIM>& getSectionBeads() {
        return *sectionBeads;
    }

    const virtual Beads<NDIM>& getSectionBeads() const {
        return *sectionBeads;
    }

    virtual Beads<NDIM>& getSectionBeads(int i) {
        return *sectionBeads;
    }

    const virtual Beads<NDIM>& getSectionBeads(int i) const {
        return *sectionBeads;
    }

    using SectionSamplerInterface::getMovingBeads;

    virtual Beads<NDIM>& getMovingBeads() {
        return *movingBeads;
    }

    const virtual Beads<NDIM>& getMovingBeads() const {
        return *movingBeads;
    }

    virtual Beads<NDIM>& getMovingBeads(int i) {
        return *movingBeads;
    }

    const virtual Beads<NDIM>& getMovingBeads(int i) const {
        return *movingBeads;
    }

    using SectionSamplerInterface::getMovingIndex;

    const virtual IArray& getMovingIndex() const {
        return *movingIndex;
    }

    const virtual IArray& getMovingIndex(int i) const {
        return *movingIndex;
    }

    const virtual SuperCell& getSuperCell() const {
        return *cell;
    }

    virtual bool isSamplingBoth() const {
        return false;
    }

    virtual AccRejEstimator* getAccRejEstimator(const std::string& name);

    virtual int getFirstSliceIndex() const {
        return sectionChooser.getFirstSliceIndex();
    }


protected:
    const int nlevel;
    const int npart;
    Beads<NDIM> *sectionBeads;
    Beads<NDIM> *movingBeads;
    IArray *movingIndex;
    IArray identityIndex;
    SectionChooser& sectionChooser;

    SuperCell* cell;
    Action *action;
    const int nrepeat;
    const int maximumMovingCount;
    CollectiveSectionMover* mover;
    AccRejEstimator* accRejEst;

};
#endif
