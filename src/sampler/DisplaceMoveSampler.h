// $Id$
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#ifndef __DisplaceMoveSampler_h_
#define __DisplaceMoveSampler_h_
class SuperCell;
class Action;
class UniformMover;
class Paths;
class ParticleChooser;
class PermutationChooser;
class Permutation;
class AccRejEstimator;
class MPIManager;

#include "algorithm/Algorithm.h"
#include <vector>
#include <cstdlib>
#include <blitz/array.h>
#include <iostream>

/** Class to perform classical displacements of the particles.
  @author Saad Khairallah, John Shumway
*/
class DisplaceMoveSampler : public Algorithm {
public:
  typedef blitz::Array<int,1> IArray;
  typedef blitz::TinyVector<double,NDIM> Vec;
  typedef blitz::Array<Vec,1> VArray;
  /// Constructor.
  DisplaceMoveSampler(const int nmoving, const int nrepeat,
    Paths&, ParticleChooser&, const UniformMover&, Action*, 
    const MPIManager* mpi);
  /// Destructor.
  virtual ~DisplaceMoveSampler();
  /// Run method, performs nrepeat samplings with probability ifreq.
  virtual void run();
  /// Get a pointer to the accept/reject statistic estimator.
  /// (You are responsible for deleting this new object.) 
  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);
  Permutation getGlobalPermutation();
  void handleBoundary(int bd1, int bd2, int sign);
protected:
  /// The number of particles to move together.
  const int nmoving;
  /// The number of times to try dispace moves.
  const int nrepeat;
  /// The number of particles.
  const int npart;
  /// The first slice this worker is moving. 
  int iFirstSlice;
  /// The last slice this worker is moving. 
  int iLastSlice;
  /// The dispalcement of the moving particles in a particular attempt.
  VArray displacement;
  /// Index of the particles to be moved.
  IArray movingIndex;
  /// Helper class to choose the particles to move.
  ParticleChooser& particleChooser; 
  /// Helper class to choose the displacement.
  const UniformMover& mover;
  /// A reference to the paths.
  Paths& paths;
  /// A reference to the simulation super cell.
  const SuperCell& cell; 
  /// The action to be evaluated during the move.
  Action *action;
  /// A pointer to the accept-reject estimator.
  AccRejEstimator* accRejEst;
  /// A pointer to the MPI manager, zero if MPI is not used.
  const MPIManager* mpi;
  /// Method to atempt a Monte Carlo move, return true if accepted.
  virtual bool tryMove();
  //  IArray * 
  const int nworker;
  blitz::Array<int,2> iworkerPerm;
 
  //Permutation localPermutation;
};
#endif
