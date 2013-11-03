#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "MPIManager.h"
#include <cstdlib>
#include <iostream>


MPIManager::MPIManager(const int nworker, const int nclone)
  : nworker(nworker), nclone(nclone)  {
#ifdef ENABLE_MPI
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
  if (rank==0 && nworker*nclone!=size) {
    std::cout << "ERROR: nworker=" << nworker
              << " and nclone=" << nclone
              << " but MPI size=" << size << std::endl;
    exit(-1);
  }
  isMainFlag=(rank==0);
  workerID=rank%nworker;
  cloneID=rank/nworker;
  workerComm=MPI::COMM_WORLD.Split(cloneID,workerID);
  int range[1][3];
  range[0][0]=0;
  range[0][1]=(nclone-1)*nworker;
  range[0][2]=nworker;
  cloneComm=MPI::COMM_WORLD.Create(
              MPI::COMM_WORLD.Get_group().Range_incl(1,range));
  isCloneMainFlag=(workerID==0);
#endif
}
