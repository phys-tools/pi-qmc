#ifndef __UniformMover_h_
#define __UniformMover_h_
#include <cstdlib>
#include <blitz/array.h>
#include <vector>

class MPIManager;
class DisplaceMoveSampler;
class DoubleDisplaceMoveSampler;

class UniformMover {
public:
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::Array<double, 1> Array;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<Vec, 1> VArray;

    UniformMover(const Vec dist, const MPIManager *mpi);
    virtual ~UniformMover();
    virtual double makeMove(VArray&, const IArray&) const;
protected:
    const MPIManager* mpi;
private:
    const Vec dist;
};
#endif
