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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "PairChooser.h"
#include "WalkingChooser.h"
#include "util/Permutation.h"
#include <iostream>

PairChooser::PairChooser(const int npart, 
  const Species &s1, const Species &s2, const int nlevel, 
  const SimulationInfo &simInfo)
  : ParticleChooser(2*npart), PermutationChooser(2*npart),
    npart(npart),
    chooser1(*new WalkingChooser(npart,s1,nlevel,simInfo)),
    chooser2(*new WalkingChooser(npart,s2,nlevel,simInfo)) {
}

PairChooser::~PairChooser() {
  delete &chooser1;
  delete &chooser2;
}

void PairChooser::chooseParticles() {
}

bool PairChooser::choosePermutation() {
  bool newPerm = chooser1.choosePermutation() &&
                 chooser2.choosePermutation();
  if (newPerm) {
    for (int i=0; i<npart; ++i) {
      index(i) = chooser1[i];
      index(i+npart) = chooser2[i];
      (*permutation)[i] = chooser1.getPermutation()[i];
      (*permutation)[i+npart] = chooser2.getPermutation()[i]+npart;
    }
  }
  return newPerm;
}

double PairChooser::getLnTranProb() const {
  return chooser1.getLnTranProb()+chooser2.getLnTranProb();
}

void PairChooser::init() {
  chooser1.init();
  chooser2.init();
}

void PairChooser::setMLSampler(const MultiLevelSampler *sampler) {
  chooser1.setMLSampler(sampler);
  chooser2.setMLSampler(sampler);
}
