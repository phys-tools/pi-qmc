#ifndef __Distance_h_
#define __Distance_h_
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

typedef blitz::TinyVector<double, NDIM> Vec;
/// Base class for one-particle distance functions.
class Distance {
public:
    virtual ~Distance() {}
    virtual double operator()(const Vec &r) {
        return 0.;
    }
};
/// Distance taken from radial separation.
class Radial: public Distance {
public:
    Radial(int idir = -1) :
            mask(1.0) {
        if (idir != -1)
            mask(idir) = 0;
    }
    virtual double operator()(const Vec &r) {
        double radius2 = 0;
        for (int i = 0; i < NDIM; ++i)
            radius2 += r(i) * r(i) * mask(i);
        return sqrt(radius2);
    }
    Vec mask;
};
/// Distance taken from Cartesian position of particle.
class Cart: public Distance {
public:
    Cart(int idim) :
            idim(idim) {
    }
    ;
    int idim;
    virtual double operator()(const Vec &r) {
        return r[idim];
    }
    ;
};

#endif
