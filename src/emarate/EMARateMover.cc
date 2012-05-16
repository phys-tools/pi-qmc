#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "EMARateMover.h"
#include "Beads.h"
#include "sampler/MultiLevelSampler.h"
#include "util/RandomNumGenerator.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "util/SuperCell.h"
#include "SimulationInfo.h"
#include "util/PeriodicGaussian.h"
#include <cmath>

EMARateMover::EMARateMover(double tau, double mass1, double mass2,
        int maxlevel, double C)
:   PermutationChooser(2),
    ParticleChooser(2),
    tau(tau),
    lambda1(0.5/mass1),
    lambda2(0.5/mass2),
    C(C) {
}

EMARateMover::~EMARateMover() {
}

double EMARateMover::calculateTransitionProbability(int nStride,
        const Beads<NDIM> &movingBeads, const Beads<NDIM> &sectionBeads,
        int nSlice, const SuperCell& cell) const {
    // Calculate 1/transition probability for move
    double oldDiagAction = 0.0;
    double oldRadAction = 0.0;
    double newDiagAction = 0.0;
    double newRadAction = 0.0;
    Vec mass1 = 0.5 / lambda1;
    Vec mass2 = 0.5 / lambda2;
    const Vec inv2Sigma21 = 0.5 * mass1 / (tau * nStride);
    const Vec inv2Sigma22 = 0.5 * mass2 / (tau * nStride);
    Vec rePrev = movingBeads(0, 0);
    Vec rhPrev = movingBeads(1, 0);
    Vec reRadPrev = rePrev;
    Vec rhRadPrev = rhPrev;
    Vec rePrevOld = sectionBeads(0, 0);
    Vec rhPrevOld = sectionBeads(1, 0);
    Vec reRadPrevOld = rePrevOld;
    Vec rhRadPrevOld = rhPrevOld;
    for (int islice = nStride; islice < nSlice; islice += nStride) {

        // Calculate action for moving beads.
        Vec re = movingBeads(0,islice);
        Vec rh = movingBeads(1,islice);

        Vec reRad = re;
        Vec rhRad = (islice == nSlice/2) ? re : rh;

        Vec delta = re - rePrev; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            newDiagAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }
        delta = rh - rhPrev; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            newDiagAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }
        delta = reRad - reRadPrev; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            newRadAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }
        delta = rhRad - rhRadPrev; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            newRadAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }
        // Calculate action for old beads.
        Vec reOld = sectionBeads(0,islice);
        Vec rhOld = sectionBeads(1,islice);

        Vec reRadOld = reOld;
        Vec rhRadOld = (islice == nSlice / 2) ? reOld : rhOld;

        delta = reOld - rePrevOld; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            oldDiagAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }
        delta = rhOld - rhPrevOld; cell.pbc(delta);
        for (int idim = 0; idim < NDIM; ++idim) {
            oldDiagAction += delta[idim] * delta[idim] * inv2Sigma22[idim];
        }
        delta = reRadOld - reRadPrevOld; cell.pbc(delta);
        for (int idim=0;idim<NDIM;++idim) {
            oldRadAction += delta[idim] * delta[idim] * inv2Sigma21[idim];
        }
        delta = rhRadOld - rhRadPrevOld; cell.pbc(delta);
        for (int idim=0;idim<NDIM;++idim) {
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
    double oldAction =
            oldDiagAction - log(1 + C * exp(-oldRadAction + oldDiagAction));
    double newAction =
            newDiagAction - log(1 + C * exp(-newRadAction + newDiagAction));
    if (C > 0.0) {
        if (log(C) - oldRadAction + oldDiagAction > 140) {
            oldAction = oldRadAction - log(C);
        }
        if (log(C) - newRadAction + newDiagAction > 140) {
            newAction =  newRadAction - log(C);
        }
    }

    return newAction - oldAction;
}

double EMARateMover::makeMove(MultiLevelSampler& sampler, const int level) {
    const Beads<NDIM>& sectionBeads = sampler.getSectionBeads();
    Beads<NDIM>& movingBeads = sampler.getMovingBeads();
    const SuperCell& cell = sampler.getSuperCell();
    const int nStride = (int)pow(2,level);
    const int nSlice = sectionBeads.getNSlice();
    double toldOverTnew = 0;
    forwardProb = 0;

    if (nStride+1 == nSlice) {
        earlierTransitions = 0.;
        isSamplingRadiating =
                chooseDiagonalOrRadiating(movingBeads, nSlice, cell);
    } else {
        moveBeads(nStride, nSlice, movingBeads, cell);

        toldOverTnew = calculateTransitionProbability(nStride, movingBeads,
                sectionBeads, nSlice, cell);
    }

    double returnValue = toldOverTnew - earlierTransitions;
    earlierTransitions = toldOverTnew;
    return returnValue;
}


double EMARateMover::calculateRadiatingProbability(
        const Beads<NDIM> & movingBeads, int nSlice, const SuperCell &cell) {
    Vec re1 = movingBeads(0, 0);
    Vec re2 = movingBeads(0, nSlice - 1);
    Vec rh1 = movingBeads(1, 0);
    Vec rh2 = movingBeads(1, nSlice - 1);
    Vec delta1 = re1 - rh1;
    cell.pbc(delta1);
    Vec delta2 = re2 - rh2;
    cell.pbc(delta2);
    int nStride = nSlice - 1;
    double inv2sigma21 = 1. / ( (lambda1 + lambda2) * tau * nStride * 2);
    double t1 = dot(delta1, delta1) * inv2sigma21;
    double t2 = dot(delta2, delta2) * inv2sigma21;
    return exp(-(t1 + t2));
}

double EMARateMover::calculateDiagonalProbability(
        const Beads<NDIM> &movingBeads, int nSlice, const SuperCell &cell) {
    Vec re1 = movingBeads(0, 0);
    Vec re2 = movingBeads(0, nSlice - 1);
    Vec rh1 = movingBeads(1, 0);
    Vec rh2 = movingBeads(1, nSlice - 1);
    Vec deltae = re2 - re1;
    cell.pbc(deltae);
    Vec deltah = rh2 - rh1;
    cell.pbc(deltah);
    int nStride = nSlice - 1;
    double inv2sigma2e = 0.25 / (lambda1 * tau * nStride);
    double inv2sigma2h = 0.25 / (lambda2 * tau * nStride);
    double te = dot(deltae, deltae) * inv2sigma2e;
    double th = dot(deltah, deltah) * inv2sigma2h;
    return exp(-(te + th));
}

bool EMARateMover::chooseDiagonalOrRadiating(
        const Beads<NDIM> & movingBeads, int nSlice,
        const SuperCell &cell) {
    double pRad = calculateRadiatingProbability(movingBeads, nSlice, cell);
    double pDiag = calculateDiagonalProbability(movingBeads, nSlice, cell);
    double randomNumber = getRandomNumber();
    return randomNumber > pDiag / (pDiag + C * pRad);
}

void EMARateMover::moveBeads(int nStride, int nSlice,
        Beads<NDIM> &movingBeads, const SuperCell &cell) {
    if (isSamplingRadiating) {
        sampleRadiating(nStride, nSlice, movingBeads);
    } else {
        sampleDiagonal(nStride, nSlice, movingBeads, cell);
    }
}

void EMARateMover::makeGaussianRandomNumbers(blitz::Array<Vec,1> & gaussRand) {
    RandomNumGenerator::makeGaussRand(gaussRand);
}

void EMARateMover::sampleRadiating(int nStride, int nSlice,
        Beads<NDIM> &movingBeads) {
    for (int iSlice = nStride; iSlice < nSlice; iSlice += 2 * nStride) {
        blitz::Array<Vec,1> gaussRand(2);
        Vec mass1 = 0.5 / lambda1;
        Vec mass2 = 0.5 / lambda2;
        Vec prev1 = (iSlice - nStride <= nSlice / 2)
                    ? movingBeads(0, iSlice - nStride)
                            : movingBeads(1, nSlice  - 1 - (iSlice - nStride));
        Vec prev2 = (iSlice - nStride < nSlice / 2)
                    ? movingBeads(0, nSlice - 1 - (iSlice - nStride))
                            : movingBeads(1, iSlice - nStride);
        Vec next1 = (iSlice + nStride <= nSlice / 2)
                    ? movingBeads(0, iSlice + nStride)
                            : movingBeads(1, nSlice -1 - (iSlice + nStride));
        Vec next2 = (iSlice + nStride < nSlice / 2)
                    ? movingBeads(0, nSlice -1 - (iSlice + nStride))
                            : movingBeads(1, iSlice + nStride);

        Vec midpoint1, midpoint2;
        // Should be rewritten to use PBC.
        if ((iSlice - nStride < nSlice / 2) && (iSlice + nStride > nSlice / 2)) {
            midpoint1 = (mass1 * prev1 + mass2 * next1) / (mass1 + mass2);
            midpoint2 = (mass1 * prev2 + mass2 * next2) / (mass1 + mass2);
        } else {
            midpoint1 = 0.5 * (prev1 + next1);
            midpoint2 = 0.5 * (prev2 + next2);
        }
        Vec sigma;
        if (iSlice - nStride < nSlice / 2) {
            if (iSlice + nStride > nSlice / 2) {
                sigma = sqrt(0.5 * (lambda1 + lambda2) * tau * nStride);
            } else {
                sigma = sqrt(lambda1 * tau * nStride);
            }
        }
        else{
            sigma = sqrt(lambda2 * tau * nStride);
        }
        makeGaussianRandomNumbers(gaussRand);
        Vec delta1 = gaussRand(0) * sigma;
        Vec delta2 = gaussRand(1) * sigma;
        if (iSlice < nSlice / 2){
            movingBeads(0, iSlice) = midpoint1 + delta1;
            movingBeads(0, nSlice - 1 - iSlice) = midpoint2 + delta2;
        } else {
            if (iSlice == nSlice / 2){
                movingBeads(0, iSlice) = midpoint1 + delta1;
                movingBeads(1, iSlice) = midpoint2 + delta2;
            } else {
                movingBeads(1, iSlice) = midpoint2 + delta2;
                movingBeads(1, nSlice - 1 - iSlice) = midpoint1 + delta1;
            }
        }
    }
}

void EMARateMover::sampleDiagonal(int nStride, int nSlice,
        Beads<NDIM> &movingBeads, const SuperCell & cell) {
    for (int iSlice = nStride; iSlice < nSlice; iSlice += 2 * nStride) {
        blitz::Array<Vec,1> gaussRand(2);
        makeGaussianRandomNumbers(gaussRand);
        for(int iMoving = 0;iMoving < 2;++iMoving){
            double lambda = (iMoving == 0) ? lambda1 : lambda2;
            double sigma = sqrt(lambda * tau * nStride);
            // Calculate the new position.
            Vec midpoint = movingBeads.delta(iMoving, iSlice + nStride, -2 * nStride);
            cell.pbc(midpoint) *= 0.5;
            midpoint += movingBeads(iMoving, iSlice - nStride);
            Vec delta = gaussRand(iMoving);
            delta *= sigma;
            movingBeads(iMoving, iSlice) = midpoint + delta;
            cell.pbc(movingBeads(iMoving, iSlice));
        }
    }
}

double EMARateMover::getRandomNumber() const {
    return RandomNumGenerator::getRand();
}

double EMARateMover::makeDelayedMove(MultiLevelSampler& sampler,
        const int level) {
    return 1.0;
}
