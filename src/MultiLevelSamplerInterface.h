#ifndef __MultiLevelSamplerInterface_h_
#define __MultiLevelSamplerInterface_h_
template <int TDIM> class Beads;
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
class SuperCell;
class SectionChooser;
#include <blitz/array.h>

class MultiLevelSamplerInterface {
public:
    typedef blitz::Array<int,1> IArray;
    virtual const Beads<NDIM>& getSectionBeads() const = 0;
    virtual const Beads<NDIM>& getMovingBeads() const = 0;
    virtual const IArray& getMovingIndex() const = 0;
    virtual int getFirstSliceIndex() const = 0;
    virtual const SuperCell& getSuperCell() const = 0;
};
#endif
