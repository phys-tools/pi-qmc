// $Id: SpinModelState.h 338 2010-11-30 18:56:16Z john.shumwayjr $
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
#ifndef __SpinModelState_h_
#define __SpinModelState_h_

#include "ModelState.h"
#include <cstdlib>
#include <blitz/array.h>

///Model state for spin up and spin down.
/// @version $Revision: 338 $
/// @author John Shumway and Jianheng Liu
class SpinModelState : public ModelState {
public:
  typedef blitz::Array<int,1> IArray;
  SpinModelState(int npart);
  virtual ~SpinModelState() {};
  virtual void write(std::ostream &os) const;
  virtual bool read(const std::string &line);
  const IArray& getSpinState() const {return spinState;}
  IArray& getModelState() {return spinState;}
  virtual int getModelCount() const {return npart+1;}
  virtual int getModelState() const;
private:
  const int npart;
  IArray spinState;
};
#endif
