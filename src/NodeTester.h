// $Id: NodeTester.h,v 1.1 2006/12/02 11:50:59 jshumwa Exp $
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifndef __NodeTester_h_
#define __NodeTester_h_

#include "Algorithm.h"
#include "LinkSummable.h"
#include <string>
template <int TDIM> class Beads;
class Paths;
class NodeModel;
class MPIManager;
class BeadFactory;

/// 
/// @version $Revision: 1.1 $
/// @author John Shumway
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
