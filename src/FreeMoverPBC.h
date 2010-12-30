// $Id$
/*  Copyright (C) 2004-2006 John B. Shumway, Jr. and Saad A. Khairallah

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
#ifndef __FreeMoverPBC_h_
#define __FreeMoverPBC_h_
#include <blitz/array.h>
#include <vector>
#include "Mover.h"
class SimulationInfo;
class PeriodicGaussian;
/** Select a trial move for beads by free particle sampling.
  * This class treats periodic boundary conditions more carefully, which
  * is necessary for very small box sizes.
  *
  * The difficulty, is that the bisection algorithm has 2**NDIM midpoints
  * when you consider periodic boundary conditions. If these are not sampled
  * carefully, some windings or permutations may be missed, raising the energy.
  *
  * @version $Revision$
  * @bug Delayed rejection move not properly coded.
  * @author John Shumway. */
class FreeMoverPBC : public Mover {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  typedef blitz::Array<double,1> Array;
  typedef blitz::Array<PeriodicGaussian*,3> PGArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::TinyVector<int,NDIM> IVec;
  /// Construct by providing lambda (@f$\lambda=1/2m@f$) timestep @f$\tau@f$.
  FreeMoverPBC(const double lambda, const int npart, const double tau);
  /// Construct by providing simInfo.
  FreeMoverPBC(const SimulationInfo&, const int maxlevel, const double pgDelta);
  /// Virtual destructor.
  virtual ~FreeMoverPBC();
  /// Move the samplers moving beads for a given level, returning
  /// the probability for the old move divided by the probability for the
  /// new move.
  virtual double makeMove(MultiLevelSampler&, const int level);  
  virtual double makeDelayedMove(MultiLevelSampler&, const int level) ;
  virtual double getForwardProb(){return forwardProb;}
private:
  /// The inverse mass, @f$\lambda=1/2m@f$.
  Array lambda;
  /// The timestep.
  double tau;
  /// Periodic gaussians.
  PGArray pg;
  /// Vector with supercell dimensions.
  Vec length; 
  /// Index of species types.
  IArray specIndex; 
  /// forward transition prob
  double forwardProb;
};
#endif
