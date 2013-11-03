#ifndef __ExchangeMover_h_
#define __ExchangeMover_h_
#include <cstdlib>
#include <blitz/array.h>
#include <vector>
#include "UniformMover.h"

class Paths;

class ExchangeMover: public UniformMover {
public:
    /// Typedefs.
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<Vec, 1> VArray;

    ExchangeMover(Paths& paths, const int& iFirstSlice, const Vec dist,
            const MPIManager *mpi);
    virtual ~ExchangeMover();
    virtual double makeMove(VArray&, const IArray&) const;
private:
    Paths& paths;
    const int iFirstSlice;
};
#endif
