#include "MultiLevelSamplerFake.h"

#include "Species.h"
#include "SimulationInfo.h"
#include "sampler/SectionSamplerInterface.h"
#include "Beads.h"
#include "util/SuperCell.h"
#include "sampler/SectionChooser.h"


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

Beads<NDIM> & MultiLevelSamplerFake::getSectionBeads() {
    return *sectionBeads;
}

Beads<NDIM> & MultiLevelSamplerFake::getMovingBeads() {
    return *movingBeads;
}

const Beads<NDIM> & MultiLevelSamplerFake::getSectionBeads() const {
    return *sectionBeads;
}

const Beads<NDIM> & MultiLevelSamplerFake::getMovingBeads() const {
    return *movingBeads;
}

const SectionSamplerInterface::IArray & MultiLevelSamplerFake::getMovingIndex() const {
    return *movingIndex;
}

Beads<NDIM> & MultiLevelSamplerFake::getSectionBeads(int i) {
    return *sectionBeads;
}

Beads<NDIM> & MultiLevelSamplerFake::getMovingBeads(int i) {
    return *movingBeads;
}

const Beads<NDIM> & MultiLevelSamplerFake::getSectionBeads(int i) const {
    return *sectionBeads;
}

const Beads<NDIM> & MultiLevelSamplerFake::getMovingBeads(int i) const {
    return *movingBeads;
}

const SectionSamplerInterface::IArray & MultiLevelSamplerFake::getMovingIndex(int i) const {
    return *movingIndex;
}


const SuperCell & MultiLevelSamplerFake::getSuperCell() const {
    return *superCell;
}


int MultiLevelSamplerFake::getFirstSliceIndex() const {
    return 0;
}
