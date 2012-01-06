#ifndef COLLECTIVESECTIONMOVER_H_
#define COLLECTIVESECTIONMOVER_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <blitz/tinyvec.h>

class CollectiveSectionMover {
public:
    typedef blitz::TinyVector<double, NDIM> Vec;

    CollectiveSectionMover(double radius, Vec amplitude);
    Vec calcShift(const Vec&) const;

    Vec getAmplitude() const;
    double getRadius() const;
    void setAmplitude(Vec amplitude);
    void setRadius(double radius);
private:
    double radius;
    Vec amplitude;
};


#endif
