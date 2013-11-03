#ifndef __NodeTester_h_
#define __NodeTester_h_

#include "algorithm/Algorithm.h"
#include "base/LinkSummable.h"
#include <string>
template <int TDIM> class Beads;
class Paths;
class NodeModel;
class MPIManager;
class BeadFactory;

class NodeTester : public Algorithm, LinkSummable {
public:
  /// Construct by providing the Paths and NodeModel to write.
  NodeTester(Paths&, NodeModel&, const std::string&, MPIManager*, 
             const BeadFactory&);
  /// Virtual destructor.
  virtual ~NodeTester() {}
  /// Write the paths.
  virtual void run();
private:
  /// The filename to write to.
  std::string filename;
  /// The Paths to use for the test.
  Paths& paths;
  /// The NodeModel to test.
  NodeModel& nodes;
  /// The MPI manager.
  MPIManager *mpi;
  /// The BeadFactory.
  const BeadFactory& beadFactory;
};
#endif
