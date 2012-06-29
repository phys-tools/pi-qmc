#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "stats/MPIManager.h"
#include "DoubleDisplaceMoveSampler.h"
#include "stats/AccRejEstimator.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "UniformMover.h"
#include "util/RandomNumGenerator.h"
#include "Paths.h"
#include "util/Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include "util/SuperCell.h"
#include <sstream>
#include <string>

DoubleDisplaceMoveSampler::DoubleDisplaceMoveSampler(int nmoving, int nrepeat,
        Paths& paths, ParticleChooser& particleChooser,
        const UniformMover& mover, Action* action, DoubleAction* doubleAction,
        const MPIManager* mpi) :
        DisplaceMoveSampler(nmoving, nrepeat, paths, particleChooser, mover,
                action, mpi), doubleAction(doubleAction), nsliceOver2(
                paths.getNSlice() / 2) {
    iFirstSlice = (paths.getLowestOwnedSlice(true));
    iLastSlice = (paths.getHighestOwnedSlice(true));
    std::cout << "In DoubleDisplaceMove :: reassign iFirstSlice = "
            << iFirstSlice << std::endl;
    std::cout << "In DoubleDisplaceMove :: reassign iLastSlice  = "
            << iLastSlice << std::endl;
    std::cout << "In DoubleDisplaceMove :: nsliceOver2 = " << nsliceOver2
            << std::endl;
}

DoubleDisplaceMoveSampler::~DoubleDisplaceMoveSampler() {
}

bool DoubleDisplaceMoveSampler::tryMove() {
    int workerID = (mpi) ? mpi->getWorkerID() : 0;
    if (workerID == 0)
        accRejEst->tryingMove(0);

    double logTranProb = mover.makeMove(displacement, movingIndex);
    if (nworker > 1) {
        handleBoundary(iFirstSlice - 1, iLastSlice + 1, +1);
        handleBoundary(iFirstSlice - 1 + nsliceOver2,
                iLastSlice + 1 + nsliceOver2, +1);
    }
    // Evaluate the change in action.
    double deltaAction =
            (action == 0) ?
                    0 :
                    action->getActionDifference(paths, displacement, nmoving,
                            movingIndex, iFirstSlice, iLastSlice);
//  //Now do the other double section
    deltaAction +=
            (action == 0) ?
                    0 :
                    action->getActionDifference(paths, displacement, nmoving,
                            movingIndex, iFirstSlice + nsliceOver2,
                            iLastSlice + nsliceOver2);
    //Now consider nodal action
    deltaAction +=
            (doubleAction == 0) ?
                    0 :
                    doubleAction->getActionDifference(paths, displacement,
                            nmoving, movingIndex, iFirstSlice, iLastSlice);

#ifdef ENABLE_MPI
    if (mpi && (mpi->getNWorker())>1) {
        double netDeltaAction=0;
        mpi->getWorkerComm().Reduce(&deltaAction,&netDeltaAction,
                1,MPI::DOUBLE,MPI::SUM,0);
        double acceptProb = exp(-netDeltaAction);
        bool acceptReject = RandomNumGenerator::getRand()>acceptProb;
        mpi->getWorkerComm().Bcast(&acceptReject,sizeof(bool),MPI::CHAR,0);
        if (acceptReject) {
            if (nworker>1) {
                handleBoundary(iFirstSlice-1, iLastSlice+1, -1);
                handleBoundary(iFirstSlice-1+nsliceOver2, iLastSlice+1+nsliceOver2, -1);
            }
            return false;
        }
    } else {
        double acceptProb=exp(-deltaAction + logTranProb);
        if (RandomNumGenerator::getRand()>acceptProb) {
            if (nworker>1) {
                handleBoundary(iFirstSlice-1, iLastSlice+1, -1);
                handleBoundary(iFirstSlice-1+nsliceOver2, iLastSlice+1+nsliceOver2, -1);
            }
            return false;
        }
    }
#else
    double acceptProb = exp(-deltaAction + logTranProb);
    if (RandomNumGenerator::getRand() > acceptProb) {
        if (nworker > 1) {
            handleBoundary(iFirstSlice - 1, iLastSlice + 1, -1);
            handleBoundary(iFirstSlice - 1 + nsliceOver2,
                    iLastSlice + 1 + nsliceOver2, -1);
        }
        return false;
    }
#endif

    if (workerID == 0)
        accRejEst->moveAccepted(0);

    // Move accepted.
    action->acceptLastMove();

    // Put moved beads in paths beads.
    for (int islice = iFirstSlice; islice <= iLastSlice; ++islice) {
        for (int imoving = 0; imoving < nmoving; ++imoving) {
            Vec &bead(paths(movingIndex(imoving), islice));
            bead += displacement(imoving);
            bead = cell.pbc(bead);

            Vec &bead2(paths(movingIndex(imoving), islice + nsliceOver2));
            bead2 += displacement(imoving);
            bead2 = cell.pbc(bead2);
        }
    }

    return true;
}
