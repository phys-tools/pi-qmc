// $Id: WorkerShifter.h,v 1.2 2006/10/18 17:08:19 jshumwa Exp $
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
#ifndef __WorkerShifter_h_
#define __WorkerShifter_h_
#include "CompositeAlgorithm.h"
class Paths;
class MPIManager;

class WorkerShifter : public CompositeAlgorithm {
public:
  WorkerShifter(const int maxShift, Paths& paths, MPIManager*);
  virtual ~WorkerShifter();
  virtual void run();
private:
  const int maxShift;
  Paths &paths;
  MPIManager* mpi;
};
#endif
