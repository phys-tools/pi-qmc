// $Id: DisplaceMoveSampler.h 22  2009-05-18 Saad Khairallah $
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
#ifndef __DisplaceMoveSampler_h_
#define __DisplaceMoveSampler_h_
template <int TDIM> class Beads;
class SuperCell;
class Action;
class DoubleAction;
class Mover;
class Paths;
class PathsChooser;
class ParticleChooser;
class PermutationChooser;
class Permutation;
class AccRejEstimator;
class BeadFactory;
class MPIManager;

#include "Algorithm.h"
#include <vector>
#include <blitz/array.h>
#include <iostream>


class DisplaceMoveSampler : public Algorithm {
public:
  typedef blitz::Array<int,1> IArray;
  DisplaceMoveSampler(const int nmoving, Paths&,  const double dist, const double freq,
		      ParticleChooser&, Mover&, Action*, const int nrepeat,
		      const BeadFactory&, const MPIManager* mpi);
  virtual ~DisplaceMoveSampler();
  virtual void run();
  Beads<NDIM>& getMovingBeads() {return *movingBeads;}
  Beads<NDIM>& getMovingBeads() const {return *movingBeads;}
  Beads<NDIM>& getPathsBeads() {return *pathsBeads;}
  const Beads<NDIM>& getPathsBeads() const {return *pathsBeads;}
  IArray&  getMovingIndex() {return *movingIndex;}
  const IArray&  getMovingIndex() const {return *movingIndex;}
  void setAction(Action*, const int level=0); /// not needed
  const SuperCell& getSuperCell() const {return cell;}
  const Paths& getPaths() const {return paths;}
  double getDist() {return dist;}
  int getNSlice() {return nslice;}
  int getNSlice() const {return nslice;}

  /// Get a pointer to the accept/reject statistic estimator.
  /// (You are responsible for deleting this new object.) 
  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);

protected:
  bool tryMove(int imovingNonPerm);

  Beads<NDIM> *pathsBeads;
  Beads<NDIM> *movingBeads;
  //const Permutation & pathsPermutation;
  ParticleChooser& particleChooser; 
  Mover& mover;
  Action *action; // DoubleAction *doubleAction;
  Paths& paths;
  AccRejEstimator* accRejEst;

  IArray *movingIndex;
  IArray identityIndex; 
  const BeadFactory& beadFactory;
  const SuperCell& cell; 
  const MPIManager* mpi;
  const int nrepeat;
  const double dist;
  const double freq;
  const int nmoving;
 
  int nslice;
  int iFirstSlice;


};
#endif
