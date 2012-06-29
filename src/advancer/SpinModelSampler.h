#ifndef __SpinModelSampler_h_
#define __SpinModelSampler_h_
class Action;
class ActionChoiceBase;
class Paths;
class AccRejEstimator;
class MPIManager;
class SpinModelState;
class Permutation;

#include "algorithm/Algorithm.h"
#include <vector>
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>

/// Class to sample different action models.

class SpinModelSampler: public Algorithm {
public:
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<Vec, 1> VArray;

    SpinModelSampler(Paths&, Action*, ActionChoiceBase*, const MPIManager* mpi);
    virtual ~SpinModelSampler();

    virtual void run();
    /// Get a pointer to the accept/reject statistic estimator.
    /// (You are responsible for deleting this new object.)
    virtual AccRejEstimator* getAccRejEstimator(const std::string& name);
protected:
    virtual bool tryMove();
    Paths& paths;
    Action *action;
    ActionChoiceBase *actionChoice;
    SpinModelState &modelState;
    const int nmodel;
    AccRejEstimator* accRejEst;
    const MPIManager* mpi;
#ifdef ENABLE_MPI
    const int nworker;
#endif
};
#endif
