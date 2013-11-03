// $Id$
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
#ifndef __WritePaths_h_
#define __WritePaths_h_

#include "Algorithm.h"
#include <string>
template <int TDIM> class Beads;
class Paths;
class MPIManager;
class BeadFactory;
class SimulationInfo;

/// Algorithm for writting paths to an ascii file.
/// Last line of the file is the same positions as the first
/// lines, just with the net permuation.
///
/// If there are clones (see MPIManager::nclone),
/// then each clone gets its own file, with its cloneID in the name. 
///
/// To handle mulitiple workers (see MPIManager::nworker)
/// writing is done in stages. Only the first worker writes.
/// The first worker writes from its lowest sample slice
/// to half its highest stored slice. Then the beads are shifted.
/// This process continues until all beads are written.
///
/// Because of the uni-directionality of shift, beads are written
/// in reverse time order.
/// 
/// @version $Revision$
/// @bug Reading and writting are done backwards because of Paths::shift
/// @author John Shumway
class WritePaths : public Algorithm {
public:
  /// Construct by providing the paths to write.
  WritePaths(Paths&, const std::string&, int dumpFreq, int maxConfigs, 
    bool writeMovie, const SimulationInfo& simInfo, MPIManager*, 
    const BeadFactory&);
  /// Virtual destructor.
  virtual ~WritePaths() {}
  /// Write the paths.
  virtual void run();
private:
  /// The filename to write to.
  std::string filename;
  /// The paths to write out.
  Paths& paths;
  /// The MPI manager.
  MPIManager *mpi;
  /// The BeadFactory.
  const BeadFactory& beadFactory;
  const int dumpFreq;
  const SimulationInfo& simInfo;
  const int maxConfigs;
  const bool writeMovie;
  std::ofstream *movieFile;
  
};
#endif

