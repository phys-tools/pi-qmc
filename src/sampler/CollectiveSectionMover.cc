#include "CollectiveSectionMover.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

#define DGESV_F77 F77_FUNC(dgesv,DGESV)
extern "C" void DGESV_F77(const int *n, const int *nrhs,
                          const double *a, const int *lda, int *ipiv,
                          double *b, const int *ldb, int *info);


CollectiveSectionMover::CollectiveSectionMover(double radius, Vec amplitude,
        Vec center, int sliceCount)
    :   radius(radius), amplitude(amplitude), center(center),
        sliceCount(sliceCount) {
}

double CollectiveSectionMover::timeEnvelope(int sliceIndex) const {
    return 1.0 - (2.0 * sliceIndex / (sliceCount - 1.0) - 1.0)
               * (2.0 * sliceIndex / (sliceCount - 1.0) - 1.0);
}

double CollectiveSectionMover::envelope(const Vec & rin,
        int sliceIndex) const {
    if (isOutsideRadius(rin)) {
        return 0.0;
    } else {
        double deltaR2 = dot(rin-center, rin-center);
        double g = 1 - deltaR2 / (radius * radius);
        g *= timeEnvelope(sliceIndex);
        return g;
    }
}

bool CollectiveSectionMover::isOutsideRadius(const Vec &rin) const {
    return dot(rin-center, rin-center) > radius * radius;
}

CollectiveSectionMover::Vec CollectiveSectionMover::calcShift(
        const Vec &rin, int sliceIndex) const {
    Vec rout = rin + envelope(rin, sliceIndex) * amplitude;
    return rout;
}

CollectiveSectionMover::Vec CollectiveSectionMover::calcInverseShift(
        const Vec &rin, int sliceIndex) const {
    std::cout << rin << std::endl;
    Vec rout = rin - envelope(rin, sliceIndex) * amplitude;
    Vec rback = calcShift(rout, sliceIndex);
    double error2 = dot(rin-rback, rin-rback)
            / (dot(rin-rout,rin-rout) + dot(rout-rback, rout-rback));
    int niter = 0;
    while (error2 > 1e-18) {
      // Solve linear equations for Newton's method.
      ++niter;
      if (niter>20) {
        std::cout << "WARNING: Too many Newton iterations "
                  << "in CollectiveSectionMover!" << std::endl;
        break;
      }
      const int N=NDIM, ONE=1;
      int info;
      IVec ipiv;
      Mat jacobian = calcJacobian(rout, sliceIndex);
      Vec delta = rin - rback;
      DGESV_F77(&N,&ONE,jacobian.data(),&N,ipiv.data(),delta.data(),&N,&info);
      rout += delta;
      rback = calcShift(rout, sliceIndex);
      error2 = dot(rin-rback, rin-rback)
              / (dot(rin-rout,rin-rout) + dot(rout-rback, rout-rback));
      }
    return rout;
}

CollectiveSectionMover::Mat CollectiveSectionMover::calcJacobian(
        const Vec &rin, int sliceIndex) const {
    Mat matrix;
    double scale;
    if (isOutsideRadius(rin)) {
        scale = 0.0;
    } else {
        scale = timeEnvelope(sliceIndex);
    }
    double invRadius2 = 1. / (radius * radius);
    for (int i=0; i<NDIM; ++i) {
        for (int j=0; j<NDIM; ++j) {
            matrix(i,j) = (i==j) ? 1.0 : 0.0;
            matrix(i,j) -= 2 * amplitude(i) * (rin(j)-center(j))
                    * invRadius2  * scale;
        }
    }
    return matrix;
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


void CollectiveSectionMover::setRadius(double radius) {
    this->radius = radius;
}
