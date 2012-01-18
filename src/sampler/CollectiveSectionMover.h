#ifndef COLLECTIVESECTIONMOVER_H_
#define COLLECTIVESECTIONMOVER_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
class SuperCell;


class CollectiveSectionMover {
public:
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::TinyVector<int, NDIM> IVec;
    typedef blitz::TinyMatrix<double,NDIM,NDIM> Mat;

    CollectiveSectionMover(double radius, Vec amplitude, Vec center,
            int level, SuperCell*);

    Vec calcShift(const Vec&, int sliceIndex) const;
    Vec calcInverseShift(const Vec&, int sliceIndex) const;
    Mat calcJacobian(const Vec&, int sliceIndex) const;

    Vec getAmplitude() const;
    double getRadius() const;
    void setAmplitude(Vec amplitude);
    void setRadius(double radius);
    int getSliceCount() const;
    Vec getCenter() const;
    void setCenter(Vec center);

private:
    inline Vec envelope(const Vec&, int sliceIndex) const;
    inline double timeEnvelope(int sliceIndex) const;
    inline bool isOutsideRadius(const Vec &rin) const;


    double radius;
    Vec amplitude;
    Vec center;
    const int sliceCount;
    SuperCell *cell;
};


#endif
