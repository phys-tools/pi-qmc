#include "MPILifecycle.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

void MPILifecycle::initialize(int argc, char **argv) {
#ifdef ENABLE_MPI
    enabled = true;
    MPI::Init(argc,argv);
    rank=MPI::COMM_WORLD.Get_rank();
#endif
}

void MPILifecycle::finalize() {
#ifdef ENABLE_MPI
    MPI::Finalize();
#endif
}

bool MPILifecycle::enabled = false;
int MPILifecycle::rank = 0;
