//$Id: SpinStatePermutationChooser.cc 338 2010-11-30 18:56:16Z john.shumwayjr $
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "SpinStatePermutationChooser.h"
#include "util/RandomNumGenerator.h"
#include "util/Permutation.h"
#include "Paths.h"
#include "MultiLevelSampler.h"
#include "Species.h"
#include "util/SuperCell.h"
#include "Beads.h"
#include "SimulationInfo.h"
#include "util/PeriodicGaussian.h"
#include "ModelState.h"
#include "SpinModelState.h"
#include <cstdlib>
#include <blitz/tinyvec-et.h>

SpinStatePermutationChooser::SpinStatePermutationChooser(const int nsize, const Species &species,
  const int nlevel, const SimulationInfo& simInfo, ModelState& modelState)
  : WalkingChooser(nsize,species,nlevel,simInfo),
    spinState(dynamic_cast<SpinModelState&>(modelState).getSpinState()) {
}

SpinStatePermutationChooser::~SpinStatePermutationChooser() {
}

void SpinStatePermutationChooser::init() {
  // Setup the table of free particle propagotor values.
  const Beads<NDIM> &sectionBeads=WalkingChooser::multiLevelSampler->getSectionBeads();
  const SuperCell &cell=multiLevelSampler->getSuperCell();
  const int nslice=sectionBeads.getNSlice(); 


  for (int ipart=0; ipart<npart; ++ipart) {
    for (int jpart=0; jpart<npart; ++jpart) {
      Vec delta=sectionBeads(ipart+ifirst,0);
      delta-=sectionBeads(jpart+ifirst,nslice-1);
      cell.pbc(delta);
      WalkingChooser::t(ipart,jpart)=1;
      for (int idim=0; idim<NDIM; ++idim) {
        WalkingChooser::t(ipart,jpart)*=(*(WalkingChooser::pg(idim)))(fabs(delta[idim]));
      }
      // Zero out the probability between different particles.
      if (spinState(ipart) != spinState(jpart)) WalkingChooser::t(ipart,jpart) = 1e-200;
    } 
  }
 
  for (int ipart=0; ipart<npart; ++ipart) {
    WalkingChooser::cump(ipart,0)=WalkingChooser::t(ipart,0);
  }
  
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int jpart=1; jpart<npart; ++jpart) {
      WalkingChooser::cump(ipart,jpart)=WalkingChooser::cump(ipart,jpart-1)+WalkingChooser::t(ipart,jpart);
      }
  }
}
