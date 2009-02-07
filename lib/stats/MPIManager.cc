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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
