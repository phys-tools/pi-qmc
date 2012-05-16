#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "EMARateAction.h"

#include <cstdlib>
#include <blitz/tinyvec-et.h>
#include "sampler/SectionSamplerInterface.h"
#include "sampler/DisplaceMoveSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "util/SuperCell.h"
#include "SimulationInfo.h"
#include "util/PeriodicGaussian.h"
#include "sampler/SectionChooser.h"

EMARateAction::EMARateAction(const SimulationInfo& simInfo,
        const Species& species1, const Species& species2, double C)
: invTau(1./simInfo.getTau()),
  species1(species1),
  species2(species2),
  index1(species1.ifirst),
  index2(species2.ifirst),
  C(C),
  nPathSlice(simInfo.getNSlice()) {

    if (species1.anMass) {
        mass1 = *species1.anMass;
    } else {
        mass1 = species1.mass;
    }
    if (species2.anMass) {
        mass2 = *species2.anMass;
    } else {
        mass2 = species2.mass;
    }

}

EMARateAction::~EMARateAction() {
}

EMARateAction::Vec EMARateAction::getMovingPosition(
        int ipart, int islice, int nMoving, const IArray& index,
        const Beads<NDIM>& sectionBeads, const Beads<NDIM>& movingBeads) {
    Vec re = sectionBeads(ipart, islice);
    for(int imoving = 0;imoving < nMoving;++imoving){
        int thisPart = index(imoving);
        if (thisPart == ipart) {
            re = movingBeads(imoving, islice);
        }
    }

    return re;
}

double EMARateAction::getActionDifference(
    const SectionSamplerInterface& sampler, const int level) {

    const Beads<NDIM>& sectionBeads(sampler.getSectionBeads());
    const Beads<NDIM>& movingBeads(sampler.getMovingBeads());

    const SuperCell& cell = sampler.getSuperCell();
    const int nStride = 1 << level;
    const int nSlice = sectionBeads.getNSlice();
    const IArray &index = sampler.getMovingIndex();
    const int nMoving = index.size();
    // Only evaluate if we are aligned with slice 0 in the middle.
    if (!isCenteredOnSliceZero(sampler, nSlice)) return 0.;

    double oldDiagAction = 0.;
    double oldRadAction = 0.;
    double newDiagAction = 0.;
    double newRadAction = 0.;

    const Vec inv2Sigma21 = 0.5 * mass1 * invTau / nStride;
    const Vec inv2Sigma22 = 0.5 * mass2 * invTau / nStride;

    Vec rePrev = getMovingPosition(
            0, 0, nMoving, index, sectionBeads,  movingBeads);
    Vec rhPrev = getMovingPosition(
            1, 0, nMoving, index, sectionBeads,  movingBeads);
    Vec reRadPrev = rePrev;
    Vec rhRadPrev = rhPrev;
    Vec rePrevOld = sectionBeads(0,0);
    Vec rhPrevOld = sectionBeads(1,0);
    Vec reRadPrevOld = rePrevOld;
    Vec rhRadPrevOld = rhPrevOld;

    for (int islice = nStride; islice < nSlice; islice += nStride) {

        // Calculate action for moving beads.
        Vec re = getMovingPosition(
                0, islice, nMoving, index, sectionBeads,  movingBeads);
        Vec rh = getMovingPosition(
                1, islice, nMoving, index, sectionBeads,  movingBeads);


        Vec reRad = re;
        Vec rhRad = (islice == nSlice / 2) ? re : rh;

        Vec delta = re - rePrev;
        cell.pbc(delta);
        for (int idim=0; idim < NDIM; ++idim) {
            newDiagAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }

        delta = rh - rhPrev;
        cell.pbc(delta);
        for (int idim=0;idim<NDIM;++idim) {
            newDiagAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }

        delta = reRad - reRadPrev;
        cell.pbc(delta);
        for (int idim=0;idim<NDIM;++idim) {
            newRadAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }

        delta=rhRad-rhRadPrev; cell.pbc(delta);
        for (int idim=0;idim<NDIM;++idim) {
            newRadAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }

        // Calculate action for old beads.
        Vec reOld = sectionBeads(0,islice);
        Vec rhOld = sectionBeads(1,islice);

        Vec reRadOld = reOld;
        Vec rhRadOld = (islice == nSlice / 2) ? reOld : rhOld;

        delta=reOld-rePrevOld; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            oldDiagAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }
        delta=rhOld-rhPrevOld; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            oldDiagAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }
        delta=reRadOld-reRadPrevOld; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            oldRadAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }
        delta=rhRadOld-rhRadPrevOld; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            oldRadAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }

        // Set the previous positions.
        rePrev = re;
        rhPrev = rh;
        reRadPrev = (islice == nSlice / 2) ? rhPrev : rePrev;
        rhRadPrev = rhPrev;
        rePrevOld = reOld;
        rhPrevOld = rhOld;
        reRadPrevOld = (islice == nSlice / 2) ? rhPrevOld : rePrevOld;
        rhRadPrevOld = rhPrevOld;

    }

    double oldAction = -log(1 + C * exp(-oldRadAction + oldDiagAction));
    double newAction = -log(1 + C * exp(-newRadAction + newDiagAction));

    if (C > 0.) {
        if (log(C) - oldRadAction + oldDiagAction > 140) {
            oldAction = oldRadAction - log(C) - oldDiagAction;
        }
        if (log(C) - newRadAction + newDiagAction > 140) {
            newAction =  newRadAction - log(C) - newDiagAction;
        }
    }

    double deltaAction = newAction - oldAction;

    return deltaAction;
}

double EMARateAction::getActionDifference(const Paths &paths, 
        const VArray &displacement, int nmoving, const IArray &movingIndex,
        int iFirstSlice, int iLastSlice) {
    return 0; // No change in action for uniform displacements of particles.
}

double EMARateAction::getTotalAction(const Paths& paths, int level) const {
    return 0;
}

void EMARateAction::getBeadAction(const Paths& paths, int ipart, int islice,
        double &u, double &utau, double &ulambda, Vec &fm, Vec &fp) const {
    u = utau = ulambda = 0; fm = 0.; fp = 0.;
}
bool EMARateAction::isCenteredOnSliceZero(
        const SectionSamplerInterface & sampler, const int nSlice) {
    const int iFirstSlice = sampler.getFirstSliceIndex();
    return iFirstSlice + nSlice / 2 == nPathSlice;
}

