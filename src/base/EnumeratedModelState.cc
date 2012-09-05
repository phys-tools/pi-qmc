#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "EnumeratedModelState.h"
#include "stats/MPIManager.h"
#include <cstdlib>

EnumeratedModelState::EnumeratedModelState(int modelCount)
  : modelCount(modelCount), modelState(0) {
}

void EnumeratedModelState::write(std::ostream &os) const {
  os << ", state " << modelState+1 << " of " << modelCount << ".";
}

bool EnumeratedModelState::read(const std::string &line) {
  std::cout << "Checking for model state..." << std::endl;
  int i = line.find("state");
  bool didRead = false;
  if (i != -1) {
    i += 6;
    int j = line.find("of");
    modelState = atoi(line.substr(i,j-i-1).c_str());
    //modelCount = atoi(line.substr(j+3).c_str());
    // Hack to handle old and new indexing of model states.
    // New convention is to start from 1 and end line with "."
    if (line[line.size()-1]=='.') --modelState;
    std::cout << "Model state " << modelState << "." << std::endl;
    didRead = true;
  }
  return didRead;
}

void EnumeratedModelState::broadcastToMPIWorkers(const MPIManager *mpi) {
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(&modelState,1,MPI::INT,0);
  }
#endif
}
