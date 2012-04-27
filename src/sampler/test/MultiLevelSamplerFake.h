#ifndef MULTILEVELSAMPLERFAKE_H_
#define MULTILEVELSAMPLERFAKE_H_

#include "sampler/SectionSamplerInterface.h"

class MultiLevelSamplerFake : public SectionSamplerInterface {
public:
    MultiLevelSamplerFake(int npart, int nmoving, int nslice);
    virtual ~MultiLevelSamplerFake();
    virtual const Beads<NDIM>& getSectionBeads() const;
    virtual const Beads<NDIM>& getMovingBeads() const;
    virtual Beads<NDIM>& getSectionBeads();
    virtual Beads<NDIM>& getMovingBeads();
    virtual const IArray& getMovingIndex() const;
    virtual const Beads<NDIM>& getSectionBeads(int) const;
    virtual const Beads<NDIM>& getMovingBeads(int) const;
    virtual Beads<NDIM>& getSectionBeads(int);
    virtual Beads<NDIM>& getMovingBeads(int);
    virtual const IArray& getMovingIndex(int) const;
    virtual bool isSamplingBoth() const {return false;}
    virtual const SuperCell& getSuperCell() const;
    virtual int getFirstSliceIndex() const;
    //virtual const SectionChooser& getSectionChooser() const;
    const int npart;
    const int nmoving;
    const int nslice;
    Beads<NDIM> *sectionBeads;
    Beads<NDIM> *movingBeads;
    IArray *movingIndex;
    SuperCell* superCell;
};


#endif
