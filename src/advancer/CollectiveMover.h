#ifndef __CollectiveMover_h_
#define __CollectiveMover_h_
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include <vector>
#include "UniformMover.h"

class Paths;

class CollectiveMover: public UniformMover {
public:
    /// Typedefs.
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::TinyVector<int, NDIM> IVec;
    typedef blitz::Array<Vec, 1> VArray;
    typedef blitz::TinyMatrix<double, NDIM, NDIM> Mat;

    CollectiveMover(Paths& paths, const int iFirstSlice, const Vec& kvec,
            const Vec& amp, const Vec &center, const IVec &random,
            const MPIManager *mpi);
    virtual ~CollectiveMover();
    virtual double makeMove(VArray&, const IArray&) const;
    void calcShift(const Vec &r) const;
    void calcJacobian(const Vec &r) const;
    void calcInverseShift(const Vec &r) const;
private:
    Paths& paths;
    const Vec kvec, amp;
    const int iFirstSlice;
    Vec center;
    IVec randomize;
    mutable Vec amplitude;
    mutable Vec phase;
    /// The value of the most recently calculated forward shift.
    mutable Vec value;
    /// The Jacobian matrix of a forward shift.
    mutable Mat jacobian;
};
#endif
