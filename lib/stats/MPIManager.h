// $Id: MPIManager.h,v 1.2 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __MPIManager_h_
#define __MPIManager_h_
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

class MPIManager {
public:
  MPIManager(const int nworker, const int nclone);
  ~MPIManager();
  bool isMain() const {return isMainFlag;}
  bool isCloneMain() const {return isCloneMainFlag;}
  int getNWorker() const {return nworker;}
  int getNClone() const {return nclone;}
  int getWorkerID() const {return workerID;}
  int getCloneID() const {return cloneID;}
#ifdef ENABLE_MPI
  const MPI::Intracomm& getWorkerComm() const {return workerComm;}
  const MPI::Intracomm& getCloneComm() const {return cloneComm;}
#endif
private:
  const int nworker;
  const int nclone;
  int workerID;
  int cloneID;
  bool isMainFlag;
  bool isCloneMainFlag;
#ifdef ENABLE_MPI
  MPI::Intracomm workerComm;
  MPI::Intracomm cloneComm;
#endif
};
#endif
