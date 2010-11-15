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
#ifndef __DoubleMLSampler_h_
#define __DoubleMLSampler_h_
template <int TDIM> class Beads;
class SuperCell;
class DoubleAction;
class Action;
class BeadFactory;
class Mover;
class Paths;
class DoubleSectionChooser;
class ParticleChooser;
class PermutationChooser;
#include "MultiLevelSampler.h"
#include "Permutation.h"
#include <vector>
#include <blitz/array.h>
#include <iostream>

/** Class for multilevel sampling of beads. 
  * @version $Revision$
  * @author John Shumway */
class DoubleMLSampler : public MultiLevelSampler {
public:
  /// Constructor.
  DoubleMLSampler(const int nmoving, Paths&, DoubleSectionChooser&,
		  ParticleChooser&, PermutationChooser&,
		  ParticleChooser&, PermutationChooser&,
		  Mover&, Action*, DoubleAction*, const bool sampleBoth, const int nrepeat,
		  const BeadFactory&, const bool delayedRejection, const double defaultFactor,
		  double newFactor);
/// Destructor.
  virtual ~DoubleMLSampler();
  /// Run the sampler.
  virtual void run();
  /// Get reference to the moving beads.
  Beads<NDIM>& getMovingBeads(const int i) {
    return i==1?*movingBeads1:*movingBeads2;}
  /// Get const reference to the moving beads.
  const Beads<NDIM>& getMovingBeads(const int i) const {
    return i==1?*movingBeads1:*movingBeads2;}
  /// Get reference to the old beads.
  Beads<NDIM>& getSectionBeads(const int i) {
    return i==1?*sectionBeads1:*sectionBeads2;}
  /// Get const reference to the moving beads.
  const Beads<NDIM>& getSectionBeads(const int i) const {
    return i==1?*sectionBeads1:*sectionBeads2;}
  /// Get a reference to the moving bead index.
  IArray&  getMovingIndex(const int i) {
    return i==1?*movingIndex1:*movingIndex2;}
  /// Get a const reference to the moving bead index.
  const IArray&  getMovingIndex(const int i) const {
    return i==1?*movingIndex1:*movingIndex2;}
  /// Represent selecting all levels by -1.
  static const int ALL_LEVELS=-1;
  /// Are both sections being sampled?
  bool isSamplingBoth() const {return samplingBoth;}
  /// Set the active section for MultiLevelSampler
  void activateSection(const int i);
  /// Get a pointer to the accept/reject statistic estimator.
  /// (You are responsible for deleting this new object.) 
  virtual AccRejEstimator* getAccRejEstimator(const std::string& name);
private:
  /// Attempt a move on the beads.
  bool tryMove(double);
  /// Reference to all beads in the section.
  Beads<NDIM> *sectionBeads1, *sectionBeads2;
  /// Reference the permutation of the in the section.
  Permutation *sectionPermutation1, *sectionPermutation2;
  /// Storage for the moving beads.
  Beads<NDIM> *movingBeads1, *movingBeads2;
  ///Storage for the rejected beads
  Beads<NDIM> *rejectedBeads1, *rejectedBeads2;
  /// Action function.
  DoubleAction *doubleAction;
  /// Index of moving particles.
  IArray *movingIndex1,*movingIndex2;
  /// More indicies for moving particles.
  IArray pMovingIndex2; 
  /// Reference to the algorithm that selected the section.
  DoubleSectionChooser& doubleSectionChooser; 
  /// Flag for determining whether to sample both sections.
  bool samplingBoth;
  /// The permutations.
  Permutation permutation1, permutation2;
  /// The algorithm for selecting the particles to move.
  ParticleChooser& particleChooser2;
  /// The algorithm for selecting the permutation.
  PermutationChooser& permutationChooser2;
  /// Number of times to repeat.
  const int nrepeat;
};
#endif
