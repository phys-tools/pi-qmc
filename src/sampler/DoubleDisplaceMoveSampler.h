#ifndef __DoubleDisplaceMoveSampler_h_
#define __DoubleDisplaceMoveSampler_h_
class DoubleAction;

#include "DisplaceMoveSampler.h"
#include <vector>
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>

/** Class to perform classical displacements of the particles.
 @author Saad Khairallah, John Shumway
 */
class DoubleDisplaceMoveSampler: public DisplaceMoveSampler {
public:
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<Vec, 1> VArray;
    /// Constructor.
    DoubleDisplaceMoveSampler(const int nmoving, const int nrepeat, Paths&,
            ParticleChooser&, const UniformMover&, Action*, DoubleAction*,
            const MPIManager* mpi);
    /// Destructor.
    virtual ~DoubleDisplaceMoveSampler();
protected:
    /// Pointer to the double action.
    DoubleAction* doubleAction;
    /// Flag for double paths.
    //const bool isDoublePaths;
    /// Half the total slices.
    const int nsliceOver2;
    /// Method to atempt a Monte Carlo move, return true if accepted.
    virtual bool tryMove();
};
#endif
