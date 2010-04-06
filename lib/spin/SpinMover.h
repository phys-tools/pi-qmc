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
#ifndef __SpinMover_h_
#define __SpinMover_h_
#include <blitz/array.h>
#include <vector>
#include "Mover.h"
#include "FreeMover.h"
class SimulationInfo;
/** Select a trial move for beads by free particle sampling.
  * @todo May want to make this a subclass of Action, with a flag
  *       for whether to return probability or allow to cancel
  *       for efficiency.
  * @version $Revision$
  * @author John Shumway. */
class SpinMover : public Mover {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::TinyVector<double,4> SVec;
  /// Construct by providing lambda (@f$\lambda=1/2m@f$) timestep @f$\tau@f$.
  SpinMover(const double lambda, const int npart, const double tau);
  /// Construct by providing simInfo.
  SpinMover(const SimulationInfo&, const int maxlevel, const double pgDelta);
  /// Virtual destructor.
  virtual ~SpinMover();
  /// Move the samplers moving beads for a given level, returning
  /// the probability for the old move divided by the probability for the
  /// new move.
  virtual double makeMove(MultiLevelSampler&, const int level);
  virtual double makeDelayedMove(MultiLevelSampler&, const int level) {return 0;} 
  virtual double getForwardProb() {return 0;}
private:
  /// A FreeMover.
  FreeMover freeMover;
  /// The inverse mass, @f$\lambda=1/2m@f$.
  Array lambda;
  /// The timestep.
  double tau;
};
#endif
