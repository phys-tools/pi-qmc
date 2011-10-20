#ifndef __SpinStatePermutationChooser_h_
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
#define __SpinStatePermutationChooser_h_

#include "PermutationChooser.h"
#include "SpeciesParticleChooser.h"
#include <cstdlib>
#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include "WalkingChooser.h"
class MultiLevelSampler;
class SimulationInfo;
class PeriodicGaussian;
class Paths;
class ModelState;
class SpinModelState;

class SpinStatePermutationChooser : public WalkingChooser {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  SpinStatePermutationChooser(const int nsize, const Species&, 
                 const int nlevel, const SimulationInfo&, ModelState& modelState);
  virtual ~SpinStatePermutationChooser();
  virtual void init();
private:
  const IArray& spinState;
};
#endif
