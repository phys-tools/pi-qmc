#ifndef __ModelSampler_h_
#define __ModelSampler_h_
class Action;
class ActionChoiceBase;
class Paths;
class AccRejEstimator;
class MPIManager;
class EnumeratedModelState;

#include "algorithm/Algorithm.h"
#include <vector>
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>

/// Class to sample different action models.
class ModelSampler: public Algorithm {
public:
    typedef blitz::Array<int, 1> IArray;
    typedef blitz::TinyVector<double, NDIM> Vec;
    typedef blitz::Array<Vec, 1> VArray;

    ModelSampler(Paths&, Action*, ActionChoiceBase*, int target,
            const MPIManager* mpi);

    virtual ~ModelSampler();

    virtual void run();
    /// Get a pointer to the accept/reject statistic estimator.
    /// (You are responsible for deleting this new object.)
    virtual AccRejEstimator* getAccRejEstimator(const std::string& name);
protected:
    /// A reference to the paths.
    Paths& paths;
    /// The action to be evaluated during the move.
    Action *action;
    /// The action to be evaluated during the move.
    ActionChoiceBase *actionChoice;
    EnumeratedModelState &modelState;
    const int nmodel;
    /// A pointer to the accept-reject estimator.
    AccRejEstimator* accRejEst;
    /// Method to atempt a Monte Carlo move, return true if accepted.
    virtual bool tryMove();
    const int target;
    /// A pointer to the MPI manager, zero if MPI is not used.
    const MPIManager* mpi;
};
#endif
