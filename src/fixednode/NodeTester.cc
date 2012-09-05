#include "config.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "NodeTester.h"
#include "NodeModel.h"
#include "base/Beads.h"
#include "base/BeadFactory.h"
#include "base/Paths.h"
#include "stats/MPIManager.h"
#include <iostream>
#include <sstream>
#include <fstream>

NodeTester::NodeTester(Paths& paths, NodeModel& nodes, 
    const std::string& filename,
    MPIManager *mpi, const BeadFactory& beadFactory)
  : filename(filename), paths(paths), nodes(nodes), 
    mpi(mpi), beadFactory(beadFactory) {
}

void NodeTester::run() {
}
