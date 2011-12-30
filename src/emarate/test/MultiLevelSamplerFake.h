#ifndef MULTILEVELSAMPLERFAKE_H_
#define MULTILEVELSAMPLERFAKE_H_

#include "MultiLevelSamplerInterface.h"

class MultiLevelSamplerFake : public MultiLevelSamplerInterface {
public:
    MultiLevelSamplerFake();
    virtual ~MultiLevelSamplerFake();
    virtual const Beads<NDIM>& getSectionBeads() const;
    virtual const Beads<NDIM>& getMovingBeads() const;
    virtual const IArray& getMovingIndex() const;
    virtual const SuperCell& getSuperCell() const;
    virtual int getFirstSliceIndex() const;
    //virtual const SectionChooser& getSectionChooser() const;
private:
    static const int npart = 2;
    static const int nmoving = 2;
    static const int nslice = 65;
    Beads<NDIM> *sectionBeads;
    Beads<NDIM> *movingBeads;
    IArray *movingIndex;
    SuperCell* superCell;
};


#endif
