// $Id: NonZeroSectionChooser.cc 338 2010-11-30 18:56:16Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr.

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
#include "NonZeroSectionChooser.h"
#include "Action.h"
#include "Paths.h"
#include "Beads.h"
#include "Permutation.h"
#include "RandomNumGenerator.h"
#include <cmath>
#include <iostream>

NonZeroSectionChooser::NonZeroSectionChooser(const int nlevel, Paths &paths, Action &action,
  const BeadFactory &beadFactory) 
  : SectionChooser(nlevel,paths.getNPart(),paths,action,beadFactory) {
}

NonZeroSectionChooser::~NonZeroSectionChooser() {
}

void NonZeroSectionChooser::run() {
  double x=RandomNumGenerator::getRand()*(1-1e-8);
  int ilo=paths->getLowestOwnedSlice(false)-1;
  if (ilo<1) ilo=1;
  int ihi=paths->getHighestSampledSlice(beads->getNSlice()-1,false);
  if (ihi+beads->getNSlice()>=paths->getNSlice()) {
    ihi = paths->getNSlice() - beads->getNSlice();
  }
  iFirstSlice=ilo+(int)((ihi+1-ilo)*x);
  if (iFirstSlice>ihi) iFirstSlice=ihi;
  // Copy coordinates from allBeads to section Beads.
  paths->getBeads(iFirstSlice,*beads);
  permutation->reset();
  // Initialize the action.
  action->initialize(*this);
  // Run the sampling algorithm.
  CompositeAlgorithm::run();
  // Copy moved coordinates from sectionBeads to allBeads.
  paths->putBeads(iFirstSlice,*beads,*permutation);
  // Refresh the buffer slices.
  paths->setBuffers();
}
