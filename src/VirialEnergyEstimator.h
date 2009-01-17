// $Id: VirialEnergyEstimator.h,v 1.13 2007/10/26 07:04:26 jshumwa Exp $
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
#ifndef __VirialEnergyEstimator_h_
#define __VirialEnergyEstimator_h_
#include "stats/ScalarEstimator.h"
#include "LinkSummable.h"
#include "Paths.h"
#include <blitz/array.h>
class Action;
class DoubleAction;
class Paths;
class SuperCell;
class SimulationInfo;
class MPIManager;
/** Virial energy estimator. 
We're going to try to subtract the following quantity,
@f[
0 = \frac{3N(L-1)}{2\tau} + \frac{1}{2\tau}\langle
\sum_{j=1}^{L-1} (R_j-R_0)\cdot\nabla_j(S_{j+1}+S_j)\rangle,
@f]
to try to reduce the varience of the energy.

See Ceperley's review article (http://link.aps.org/abstract/RMP/v67/p279)
or Herman, Bruskin, and Berne (J. Chem. Phys. 79, 5150-5155, 1983).
  
@f[ E_V = \left\langle \frac{dU_i}{d\tau} -
                       \frac{1}{2} F_i\Delta_i \right\rangle @f].

@bug In development and testing.
@version $Revision: 1.13 $
@author John Shumway */
class VirialEnergyEstimator : public ScalarEstimator, public LinkSummable {
public:
  /// Number of dimensions.
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  typedef blitz::Array<bool,1> BArray;
  /// Constructor.
  VirialEnergyEstimator(const SimulationInfo&, const Action*, 
    const DoubleAction*, const int nwindow, MPIManager *mpi,
    const std::string& unitName, double scale, double shift);
  /// Virtual destructor.
  virtual ~VirialEnergyEstimator() {}
  /// Initialize the calculation.
  virtual void initCalc(const int nslice, const int firstSlice);
  /// Add contribution from a link.
  virtual void handleLink(const Vec& start, const Vec& end,
                          const int ipart, const int islice, const Paths&);
  /// Finalize the calculation.
  virtual void endCalc(const int nslice);
  /// Get value of virial energy estimator.
  virtual double calcValue() {return value=etot/enorm;}
  /// Clear value of virial energy estimator.
  virtual void reset() {etot=enorm=0;}
  /// Evaluate for Paths configuration.
  virtual void evaluate(const Paths& paths) {paths.sumOverLinks(*this);}
private:
  /// The timestep.
  const double tau;
  /// Virial energy estimate.
  double energy;
  /// The accumulated energy.
  double etot;
  /// The normalization of the accumualted energy.
  double enorm;
  /// The number of particles.
  const int npart;
  /// The action.
  const Action* action;
  /// The double action.
  const DoubleAction* doubleAction;
  /// The staring position of the path (needed for PBC).
  VArray r0;
  /// The average position of the path.
  VArray rav;
  /// The average force on the path.
  VArray fav;
  /// The supercell.
  SuperCell& cell;
  /// The window size.
  const int nwindow;
  /// The total number of slices.
  int nslice;
  /// The first slice.
  int firstSlice;
  /// Flag for static particles.
  BArray isStatic;
  /// The MPI Manager.
  MPIManager *mpi; 
};

#endif
