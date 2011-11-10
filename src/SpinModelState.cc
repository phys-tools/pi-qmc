// $Id: SpinModelState.cc 338 2010-11-30 18:56:16Z john.shumwayjr $
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "SpinModelState.h"
#include <cstdlib>
#include "stats/MPIManager.h"

SpinModelState::SpinModelState(int npart)
  : npart(npart), spinState(npart) {
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
