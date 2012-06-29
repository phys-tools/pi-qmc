#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "stats/MPIManager.h"
#include "ExchangeMover.h"
#include "Paths.h"
#include "Beads.h"
#include "util/RandomNumGenerator.h"
#include "util/SuperCell.h"
#include "SimulationInfo.h"
#include <cmath>
#include <sstream>
#include <string>

ExchangeMover::ExchangeMover(Paths& paths, const int& iFirstSlice,
        const Vec dist, const MPIManager* mpi) :
        UniformMover(dist, mpi), paths(paths), iFirstSlice(iFirstSlice) {
}

ExchangeMover::~ExchangeMover() {
}

double ExchangeMover::makeMove(VArray& displacement,
        const IArray& movingIndex) const {
    const Vec& bead1 = paths(movingIndex(0), iFirstSlice);
    const Vec& bead2 = paths(movingIndex(1), iFirstSlice);
    displacement(0) = bead2 - bead1;
    displacement(1) = bead1 - bead2;
#ifdef ENABLE_MPI
    if (mpi && (mpi->getNWorker())>1) {
        mpi->getWorkerComm().Bcast(displacement.data(),2*NDIM,MPI::DOUBLE,0);
    }
#endif
    return 0;
}

