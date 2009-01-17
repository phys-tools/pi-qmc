// $Id: DampedFreeTensorMover.h,v 1.4 2006/10/18 17:08:18 jshumwa Exp $
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
#ifndef __DampedFreeTensorMover_h_
#define __DampedFreeTensorMover_h_
#include "Mover.h"
class SimulationInfo;
#include <blitz/array.h>
/** Select a trial move for beads by free particle sampling.
  * @todo May want to make this a subclass of Action, with a flag
  *       for whether to return probability or allow to cancel
  *       for efficiency.
  * @bug  Needs to calculate the transiton probability correctly when
  *       periodic images are significant.  Currently the routine
  *       neglects this.
  * @bug  Hard-coded for Si/Ge electron and hole.
  * @version $Revision: 1.4 $
  * @author John Shumway. */
class DampedFreeTensorMover : public Mover {
public:
  typedef blitz::TinyVector<double,NDIM> Vec;
  /// Construct by providing lambda (@f$\lambda=1/2m@f$) timestep @f$\tau@f$.
  DampedFreeTensorMover(const SimulationInfo &simInfo,
                        const int saturationLevel);
  /// Virtual destructor.
  virtual ~DampedFreeTensorMover() {}
  /// Move the samplers moving beads for a given level, returning
  /// the probability for the old move divided by the probability for the
  /// new move.
  virtual double makeMove(MultiLevelSampler&, const int level);
private:
  /// The inverse mass, @f$\lambda=1/2m@f$.
  blitz::Array<blitz::TinyVector<double,NDIM>,1> lambda;
  /// The timestep.
  double tau;
  /// The level at which lambda saturation occurs.
  int saturationLevel;
};
#endif
