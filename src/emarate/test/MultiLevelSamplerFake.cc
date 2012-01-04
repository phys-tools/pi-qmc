#include "MultiLevelSamplerFake.h"

#include "Species.h"
#include "SimulationInfo.h"
#include "MultiLevelSamplerInterface.h"
#include "Beads.h"
#include "SuperCell.h"
#include "SectionChooser.h"


MultiLevelSamplerFake::MultiLevelSamplerFake(int npart, int nmoving, int nslice)
:   npart(npart), nmoving(nmoving), nslice(nslice),
    sectionBeads(new Beads<NDIM>(npart, nslice)),
    movingBeads(new Beads<NDIM>(nmoving, nslice)),
    movingIndex(new IArray(nmoving)),
    superCell(new SuperCell(blitz::TinyVector<double, NDIM>(1.0,0))) {
}

MultiLevelSamplerFake::~MultiLevelSamplerFake() {
    delete movingBeads;
    delete sectionBeads;
    delete movingIndex;
    delete superCell;
}

const Beads<NDIM> & MultiLevelSamplerFake::getSectionBeads() const {
    return *sectionBeads;
}



const Beads<NDIM> & MultiLevelSamplerFake::getMovingBeads() const {
    return *movingBeads;
}



const MultiLevelSamplerInterface::IArray & MultiLevelSamplerFake::getMovingIndex() const {
    return *movingIndex;
}

const SuperCell & MultiLevelSamplerFake::getSuperCell() const {
    return *superCell;
}


int MultiLevelSamplerFake::getFirstSliceIndex() const {
    return 0;
}
