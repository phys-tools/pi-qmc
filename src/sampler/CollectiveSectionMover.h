#ifndef COLLECTIVESECTIONMOVER_H_
#define COLLECTIVESECTIONMOVER_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include "SectionSamplerInterface.h"
class SuperCell;
class CollectiveSectionSampler;

class CollectiveSectionMover {
public:
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::TinyVector<int, NDIM> IVec;
    typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;

    CollectiveSectionMover(SuperCell* cell);
    ~CollectiveSectionMover();
    double makeMove(CollectiveSectionSampler&);
    Vec calcShift(const Vec&, int sliceIndex) const;
    Vec calcInverseShift(const Vec&, int sliceIndex) const;
    Mat calcJacobian(const Vec&, int sliceIndex) const;
    double calcJacobianDet(const Mat&);

    void setAmplitude(Vec amplitude);
    Vec getAmplitude() const;

    void setRadius(double radius);
    void setCenter(Vec center);

    double getRadius() const;
    void setSliceCount(int);

    int getSliceCount() const;
    Vec getCenter() const;

private:
    inline Vec calcDisplacement(const Vec&, int sliceIndex) const;
    inline double timeEnvelope(int sliceIndex) const;
    inline bool isOutsideRadius(const Vec &rin) const;

    double radius;
    Vec amplitude;
    Vec center;
    int sliceCount;
    SuperCell *cell;
};


#endif
