// $Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr.

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
#ifndef __PairChooser_h_
#define __PairChooser_h_

#include "ParticleChooser.h"
#include "PermutationChooser.h"
class MultiLevelSampler;
class Species;
class SimulationInfo;
class WalkingChooser;
/*** Choose permutations for two species.
  * @version $Revision$
  * @author John Shumway */
class PairChooser : public ParticleChooser, public PermutationChooser {
public:
  /// Constructor.
  PairChooser(const int npart, const Species&, const Species&, 
              const int nlevel, const SimulationInfo&);
  /// Virtual destructor.
  virtual ~PairChooser();
  /// Choose a new set of particles.
  virtual void chooseParticles();
  /// Choose a new permutation.
  virtual bool choosePermutation();
  /// Initialize.
  virtual void init();
  /// Get log of the transition probability.
  virtual double getLnTranProb() const;
  /// Set the MultiLevelSampler (needed by init method).
  virtual void setMLSampler(const MultiLevelSampler*);
private:
  int npart;
  WalkingChooser &chooser1, &chooser2;
};
#endif
