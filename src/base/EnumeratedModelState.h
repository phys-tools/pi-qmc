// $Id: EnumeratedModelState.h 338 2010-11-30 18:56:16Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr and Jianheng Liu.

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
#ifndef __EnumeratedModelState_h_
#define __EnumeratedModelState_h_

#include "ModelState.h"

///Base class for model state.
/// @version $Revision: 338 $
/// @author John Shumway and Jianheng Liu
class EnumeratedModelState : public ModelState {
public:
  EnumeratedModelState(int modelCount);
  virtual ~EnumeratedModelState() {};
  virtual void write(std::ostream &os) const;
  virtual bool read(const std::string &line);
  void setModelState(int i) {modelState = i;}
  virtual int getModelCount() const {return modelCount;}
  virtual int getModelState() const {return modelState;}
  virtual void broadcastToMPIWorkers(const MPIManager *mpi);
private:
  const int modelCount;
  int modelState;
};
#endif
