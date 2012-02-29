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

    CollectiveSectionMover(double radius, Vec amplitude, int npart, Vec min,
	                                           Vec max, SuperCell* cell);
    ~CollectiveSectionMover();
    double makeMove(CollectiveSectionSampler& sampler, int ilevel);
    Vec calcShift(const Vec&, int sliceIndex) const;
    Vec calcInverseShift(const Vec&, int sliceIndex) const;
    Mat calcJacobian(const Vec&, int sliceIndex) const;
    double calcJacobianDet(const Mat&);
/*
    Vec getAmplitude() const;
    double getRadius() const;
    void setAmplitude(Vec amplitude);
    void setRadius(double radius);
    int getSliceCount() const;
    Vec getCenter() const;
    void setCenter(Vec center);
*/
private:
    inline Vec calcDisplacement(const Vec&, int sliceIndex) const;
    inline double timeEnvelope(int sliceIndex) const;
    inline bool isOutsideRadius(const Vec &rin) const;

    double radius;
    Vec amplitude, amp;
    /// Center of the cylinder.
    Vec center;
    /// Boundary of for the center.
    Vec min, max;
    int sliceCount;
    const int npart;
    SuperCell *cell;
};


#endif
