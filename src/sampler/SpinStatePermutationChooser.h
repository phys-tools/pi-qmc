#ifndef __SpinStatePermutationChooser_h_
#define __SpinStatePermutationChooser_h_

#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include "WalkingChooser.h"
class SimulationInfo;
class ModelState;
class SpinModelState;
class Species;

class SpinStatePermutationChooser: public WalkingChooser {
public:
    typedef blitz::Array<int, 1> IArray;
    SpinStatePermutationChooser(const int nsize, const Species&,
            const int nlevel, const SimulationInfo&, ModelState& modelState);
    virtual ~SpinStatePermutationChooser();
    virtual void init();
private:
    const IArray& spinState;
};
#endif
