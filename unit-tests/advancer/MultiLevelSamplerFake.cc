#include "MultiLevelSamplerFake.h"
#include "advancer/SectionChooser.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Species.h"
#include "base/SimulationInfo.h"
#include "base/Beads.h"
#include "util/SuperCell.h"

MultiLevelSamplerFake::MultiLevelSamplerFake(int npart, int nmoving, int nslice)
:   npart(npart), nmoving(nmoving), nslice(nslice),
    firstSliceIndex(0),
    sectionBeads(new Beads<NDIM>(npart, nslice)),
    movingBeads(new Beads<NDIM>(nmoving, nslice)),
    movingIndex(new IArray(nmoving)) {

    blitz::TinyVector<double, NDIM> length = 10.0;
    superCell = new SuperCell(length);

    for (int ipart = 0; ipart < nmoving; ++ ipart) {
        (*movingIndex)(ipart) = ipart;
    }
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

const SectionSamplerInterface::IArray&
MultiLevelSamplerFake::getMovingIndex() const {
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

const SectionSamplerInterface::IArray&
MultiLevelSamplerFake::getMovingIndex(int i) const {
    return *movingIndex;
}


const SuperCell & MultiLevelSamplerFake::getSuperCell() const {
    return *superCell;
}


int MultiLevelSamplerFake::getFirstSliceIndex() const {
    return firstSliceIndex;
}

void MultiLevelSamplerFake::copySectionBeadsToMovingBeads() {
    for (int islice = 0;islice < nslice;++islice){
        for (int imoving = 0;imoving < nmoving;++imoving){
            int ipart = (*movingIndex)(imoving);
            (*movingBeads)(imoving, islice) = (*sectionBeads)(ipart, islice);
        }
    }
}
