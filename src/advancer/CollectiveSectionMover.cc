#include "CollectiveSectionMover.h"
#include "SectionSamplerInterface.h"
#include "CollectiveSectionSampler.h"
#include "base/Beads.h"
#include "util/RandomNumGenerator.h"
#include "util/SuperCell.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <blitz/tinyvec-et.h>

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
extern "C" void DGESV_F77(const int *n, const int *nrhs, const double *a,
        const int *lda, int *ipiv, double *b, const int *ldb, int *info);

CollectiveSectionMover::CollectiveSectionMover(SuperCell* cell) :
        cell(cell) {
    center = 0.0;
    radius = 1.0;
    sliceCount = 3;
    amplitude = 0.0;
}

CollectiveSectionMover::~CollectiveSectionMover() {
}

double CollectiveSectionMover::makeMove(CollectiveSectionSampler& sampler) {
    const Beads<NDIM>& sectionBeads = sampler.getSectionBeads();
    Beads<NDIM>& movingBeads = sampler.getMovingBeads();
    int npart = movingBeads.getNPart();
    sliceCount = sectionBeads.getNSlice();
    double tranProb = 1.;
    cell->pbc(center);
    bool forward = (RandomNumGenerator::getRand() > 0.5);
    double jacobian = 0.;
    for (int islice = 1; islice < sliceCount - 1; ++islice) {
        for (int ipart = 0; ipart < npart; ++ipart) {
            jacobian = calcJacobianDet(
                    calcJacobian(movingBeads(ipart, islice), islice));
            if (forward) {
                movingBeads(ipart, islice) = calcShift(
                        movingBeads(ipart, islice), islice);
                tranProb *= jacobian;
            } else {
                movingBeads(ipart, islice) = calcInverseShift(
                        movingBeads(ipart, islice), islice);
                tranProb /= jacobian;
            }
        }
    }
    return log(tranProb);
}

double CollectiveSectionMover::timeEnvelope(int sliceIndex) const {
    return 1.0
            - (2.0 * sliceIndex / (sliceCount - 1.0) - 1.0)
                    * (2.0 * sliceIndex / (sliceCount - 1.0) - 1.0);
}

CollectiveSectionMover::Vec CollectiveSectionMover::calcDisplacement(
        const Vec & rin, int sliceIndex) const {
    if (isOutsideRadius(rin)) {
        return 0.0;
    } else {
        Vec delta = rin - center;
        cell->pbc(delta);
        double deltaR2 = dot(delta, delta);
        double g = 1 - deltaR2 / (radius * radius);
        g *= timeEnvelope(sliceIndex);
        return g * amplitude;
    }
}

bool CollectiveSectionMover::isOutsideRadius(const Vec &rin) const {
    Vec delta = rin - center;
    cell->pbc(delta);
    return dot(delta, delta) > radius * radius;
}

CollectiveSectionMover::Vec CollectiveSectionMover::calcShift(const Vec &rin,
        int sliceIndex) const {
    Vec rout = rin + calcDisplacement(rin, sliceIndex);
    cell->pbc(rout);
    return rout;
}

CollectiveSectionMover::Vec CollectiveSectionMover::calcInverseShift(
        const Vec &rin, int sliceIndex) const {
    if (isOutsideRadius(rin))
        return rin;
    /// Calculate the inverse shift by Newton iteration.
    Vec rout = rin - calcDisplacement(rin, sliceIndex);
    Vec rback = calcShift(rout, sliceIndex);
    Vec deltaInBack = rin - rback;
    cell->pbc(deltaInBack);
    Vec deltaInOut = rin - rout;
    cell->pbc(deltaInOut);
    Vec deltaOutBack = rout - rback;
    double error2 = dot(deltaInBack, deltaInBack)
            / (dot(deltaInOut, deltaInOut) + dot(deltaOutBack, deltaOutBack));
    int niter = 0;
    while (error2 > 1e-25) {
        // Solve linear equations for Newton's method.
        ++niter;
        if (niter > 25) {
            std::cout << "WARNING: Too many Newton iterations "
                    << "in CollectiveSectionMover!" << std::endl;
            break;
        }
        const int N = NDIM, ONE = 1;
        int info;
        IVec ipiv;
        Mat jacobian = calcJacobian(rout, sliceIndex);
        Vec delta = rin - rback;
        cell->pbc(delta);
        DGESV_F77(&N, &ONE, jacobian.data(), &N, ipiv.data(), delta.data(), &N,
                &info);
        rout += delta;
        cell->pbc(rout);
        rback = calcShift(rout, sliceIndex);
        Vec deltaInBack = rin - rback;
        cell->pbc(deltaInBack);
        Vec deltaInOut = rin - rout;
        cell->pbc(deltaInOut);
        Vec deltaOutBack = rout - rback;
        error2 =
                dot(deltaInBack, deltaInBack)
                        / (dot(deltaInOut, deltaInOut)
                                + dot(deltaOutBack, deltaOutBack));
    }
    cell->pbc(rout);
    return rout;
}

CollectiveSectionMover::Mat CollectiveSectionMover::calcJacobian(const Vec &rin,
        int sliceIndex) const {
    Mat matrix;
    double scale;
    if (isOutsideRadius(rin)) {
        scale = 0.0;
    } else {
        scale = timeEnvelope(sliceIndex);
    }
    double invRadius2 = 1. / (radius * radius);
    Vec delta = rin - center;
    cell->pbc(delta);
    for (int i = 0; i < NDIM; ++i) {
        for (int j = 0; j < NDIM; ++j) {
            matrix(i, j) = (i == j) ? 1.0 : 0.0;
            matrix(i, j) -= 2 * amplitude(i) * delta(j) * invRadius2 * scale;
        }
    }
    return matrix;
}

double CollectiveSectionMover::calcJacobianDet(
        const CollectiveSectionMover::Mat& jacobian) {
#if NDIM==1
    double jacob = jacobian(0,0);
#elif NDIM==2
    double jacob = jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0);
#elif NDIM==3
    double jacob =
            jacobian(0, 0)
                    * (jacobian(1, 1) * jacobian(2, 2)
                            - jacobian(1, 2) * jacobian(2, 1))
                    - jacobian(0, 1)
                            * (jacobian(1, 2) * jacobian(2, 0)
                                    - jacobian(1, 0) * jacobian(2, 2))
                    + jacobian(0, 2)
                            * (jacobian(1, 0) * jacobian(2, 1)
                                    - jacobian(1, 1) * jacobian(2, 0));
#endif
    return jacob;
}

CollectiveSectionMover::Vec CollectiveSectionMover::getAmplitude() const {
    return amplitude;
}

double CollectiveSectionMover::getRadius() const {
    return radius;
}

void CollectiveSectionMover::setAmplitude(Vec amplitude) {
    this->amplitude = amplitude;
}

int CollectiveSectionMover::getSliceCount() const {
    return sliceCount;
}

void CollectiveSectionMover::setSliceCount(int count) {
    this->sliceCount = count;
}

void CollectiveSectionMover::setRadius(double radius) {
    this->radius = radius;
}

CollectiveSectionMover::Vec CollectiveSectionMover::getCenter() const {
    return center;
}

void CollectiveSectionMover::setCenter(Vec center) {
    this->center = center;
}

