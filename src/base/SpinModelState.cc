#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "SpinModelState.h"
#include <cstdlib>
#include "stats/MPIManager.h"

SpinModelState::SpinModelState(int npart, int initial)
  : npart(npart), spinState(npart) {
  if (initial == 1)
    for (int i=0; i<npart; ++i) spinState(i) = (i%2 == 0) ? 0 : 1;
  else
    for (int i=0; i<npart; ++i) spinState(i) = (i<npart/2) ? 0 : 1;
}

void SpinModelState::write(std::ostream &os) const {
  // Writes ", state: 0 1 0 0 1 1 0" to stream os.
  os << ", state:";
  for (int i=0; i<npart; ++i) os << " " << spinState(i);
}

bool SpinModelState::read(const std::string &line) {
  // Reads "..., state: 0 1 0 0 1 1 0" from string line.
  std::cout << "Checking for spin model state..." << std::endl;
  int i = line.find("state");
  bool didRead = false;
  if (i != -1) {
    int offset = i+7;
    for (int j=0; j<npart; ++j) {
      spinState(j) = int(line[offset+2*j]-'0');
    }
    didRead = true;
  }
  if (didRead == true)
    std::cout << "Found spin model state: " <<spinState<<std::endl;
  return didRead;
}

int SpinModelState::getModelState() const {
  int sztot =  blitz::sum(spinState);
  return sztot;
}

void SpinModelState::broadcastToMPIWorkers(const MPIManager *mpi) {
#ifdef ENABLE_MPI
  if (mpi)
    mpi->getWorkerComm().Bcast(spinState.data(),
                               npart+1,MPI::INT,0);
#endif
}
