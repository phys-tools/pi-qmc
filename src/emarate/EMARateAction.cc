#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "EMARateAction.h"
#include "advancer/DisplaceMoveSampler.h"
#include "advancer/SectionSamplerInterface.h"
#include "action/coulomb/CoulombLinkAction.h"
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "util/SuperCell.h"
#include "util/PeriodicGaussian.h"
#include <cstdlib>
#include <blitz/tinyvec-et.h>

EMARateAction::EMARateAction(const SimulationInfo& simInfo,
        const Species* species1, const Species* species2, double C)
:   invTau(1./simInfo.getTau()),
    species1(species1),
    species2(species2),
    index1(species1->ifirst),
    index2(species2->ifirst),
    C(C),
    nPathSlice(simInfo.getNSlice()),
    hasCoulomb(false),
    coulomb(0) {

    if (species1->anMass) {
        mass1 = *(species1->anMass);
    } else {
        mass1 = species1->mass;
    }
    if (species2->anMass) {
        mass2 = *(species2->anMass);
    } else {
        mass2 = species2->mass;
    }

}

EMARateAction::~EMARateAction() {
    delete coulomb;
}

void EMARateAction::includeCoulombContribution(double epsilon, int norder) {
    hasCoulomb = true;
    double q1q2 = -1.0;
    double mu = 1.0 / (1.0 / species1->mass + 1.0 / species2->mass);
    double deltaTau = 1.0 / invTau;
    coulomb = new CoulombLinkAction(q1q2, epsilon, mu, deltaTau, norder);
}

EMARateAction::Vec EMARateAction::getMovingPosition(
        int ipart, int islice, int nMoving, const IArray& index,
        const Beads<NDIM>& sectionBeads, const Beads<NDIM>& movingBeads) {
    Vec position = sectionBeads(ipart, islice);
    for (int imoving = 0;imoving < nMoving; ++imoving){
        int thisPart = index(imoving);
        if (thisPart == ipart) {
            position = movingBeads(imoving, islice);
        }
    }
    return position;
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

    double oldDiagAction = 0.0;
    double oldRadAction = 0.0;
    double newDiagAction = 0.0;
    double newRadAction = 0.0;

    const Vec inv2Sigma21 = 0.5 * mass1 * invTau / nStride;
    const Vec inv2Sigma22 = 0.5 * mass2 * invTau / nStride;

    Vec rePrev = getMovingPosition(
            index1, 0, nMoving, index, sectionBeads,  movingBeads);
    Vec rhPrev = getMovingPosition(
            index2, 0, nMoving, index, sectionBeads,  movingBeads);
    Vec reRadPrev = rePrev;
    Vec rhRadPrev = rhPrev;
    Vec rePrevOld = sectionBeads(index1,0);
    Vec rhPrevOld = sectionBeads(index2,0);
    Vec reRadPrevOld = rePrevOld;
    Vec rhRadPrevOld = rhPrevOld;

    for (int islice = nStride; islice < nSlice; islice += nStride) {

        // Calculate action for moving beads.
        Vec re = getMovingPosition(
                index1, islice, nMoving, index, sectionBeads,  movingBeads);
        Vec rh = getMovingPosition(
                index2, islice, nMoving, index, sectionBeads,  movingBeads);


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

        if (hasCoulomb && islice == nSlice/2 && nStride == 1) {
            Vec deltaPrev = reRadPrev - rhRadPrev;
            Vec deltaRad = reRad - rhRad;
            Vec deltaDiag = re - rh;
            newDiagAction += coulomb->getValue(deltaDiag, deltaPrev);
            newRadAction += coulomb->getValue(deltaRad, deltaPrev);
        }

        if (hasCoulomb && islice == nSlice/2 + nStride && nStride == 1) {
            Vec deltaRad = reRadPrev - rhRadPrev;
            Vec deltaDiag = rePrev - rhRad;
            Vec delta = re - rh;
            newDiagAction += coulomb->getValue(delta, deltaDiag);
            newRadAction += coulomb->getValue(delta, deltaRad);
        }

        // Calculate action for old beads.
        Vec reOld = sectionBeads(index1,islice);
        Vec rhOld = sectionBeads(index2,islice);

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

        if (hasCoulomb && islice == nSlice/2 && nStride == 1) {
            Vec deltaPrev = reRadPrevOld - rhRadPrevOld;
            Vec deltaRad = reRadOld - rhRadOld;
            Vec deltaDiag = reOld - rhOld;
            oldDiagAction += coulomb->getValue(deltaDiag, deltaPrev);
            oldRadAction += coulomb->getValue(deltaRad, deltaPrev);
        }

        if (hasCoulomb && islice == nSlice/2 + nStride && nStride == 1) {
            Vec deltaRad = reRadPrevOld - rhRadPrevOld;
            Vec deltaDiag = rePrevOld - rhPrevOld;
            Vec delta = reOld - rhOld;
            oldDiagAction += coulomb->getValue(delta, deltaDiag);
            oldRadAction += coulomb->getValue(delta, deltaRad);
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

    if (C > 0.0) {
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
