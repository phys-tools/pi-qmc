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
#ifndef __MultiLevelSampler_h_
#define __MultiLevelSampler_h_
template <int TDIM> class Beads;
class SuperCell;
class Action;
class Mover;
class Paths;
class SectionChooser;
class ParticleChooser;
class PermutationChooser;
class Permutation;
class AccRejEstimator;
class BeadFactory;
#include "Algorithm.h"
#include <vector>
#include <blitz/array.h>
#include <iostream>

/** Class for multilevel sampling of beads. 
  * @version $Revision$
  * @author John Shumway */
class MultiLevelSampler : public Algorithm {
public:
  /// Typedefs.
  typedef blitz::Array<int,1> IArray;
  /// Constructor.
  MultiLevelSampler(const int nmoving, Paths&, SectionChooser&,
    ParticleChooser&, PermutationChooser&, Mover&, Action*, const int nrepeat,
    const BeadFactory&);
  /// Destructor.
  virtual ~MultiLevelSampler();
  /// Run the sampler.
  virtual void run();
  /// Get the number of levels in the sampling.
  int getNLevel() const {return nlevel;}
  /// Get reference to the moving beads.
  Beads<NDIM>& getMovingBeads() {return *movingBeads;}
  /// Get const reference to the moving beads.
  const Beads<NDIM>& getMovingBeads() const {return *movingBeads;}
  /// Get reference to the old beads.
  Beads<NDIM>& getSectionBeads() {return *sectionBeads;}
  /// Get const reference to the moving beads.
  const Beads<NDIM>& getSectionBeads() const {return *sectionBeads;}
  /// Get a reference to the moving bead index.
  IArray&  getMovingIndex() {return *movingIndex;}
  /// Get a const reference to the moving bead index.
  const IArray&  getMovingIndex() const {return *movingIndex;}
  /// Represent selecting all levels by -1.
  static const int ALL_LEVELS=-1;
  /// Set the action function for a level, or default to level=ALL_LEVELS.
  void setAction(Action*, const int level=ALL_LEVELS);
  /// Get a const reference to the SuperCell. 
  const SuperCell& getSuperCell() const {return cell;}
  /// Get a constant reference to the paths.
  const Paths& getPaths() const {return paths;}
  /// Get a pointer to the accept/reject statistic estimator.
  /// (You are responsible for deleting this new object.) 
  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);
protected:
  /// Attempt a move on the beads.
  bool tryMove(double initialTranProb);
  /// Number of levels.
  const int nlevel;
  /// Number of moving particles.
  const int nmoving;
  /// Reference to all beads in the section.
  Beads<NDIM> *sectionBeads;
  /// Reference the permutation of the in the section.
  Permutation *sectionPermutation;
  /// Storage for the moving beads.
  Beads<NDIM> *movingBeads;
  /// Algorithm to select trial move for the beads.
  Mover& mover;
  /// Supercell.
  const SuperCell& cell;
  /// Action function.
  Action *action;
  /// Index of moving particles.
  IArray *movingIndex;
  /// More indicies for moving particles.
  IArray identityIndex, pMovingIndex;
  /// The algorithm for selecting the particles to move.
  ParticleChooser& particleChooser; 
  /// The algorithm for selecting the permutation.
  PermutationChooser& permutationChooser; 
  /// Reference to the algorithm that selected the section.
  SectionChooser& sectionChooser; 
  /// Reference to the paths (only needed for fixed node reference point).
  const Paths& paths;
  /// AccRejEstimator.
  AccRejEstimator* accRejEst;
  /// Number of times to repeat.
  const int nrepeat;
};
#endif
