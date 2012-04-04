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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "RandomNumGenerator.h"
#include "DoubleSectionChooser.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "Paths.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "Permutation.h"
#include <cmath>
#include <iostream>
#include <gsl/gsl_qrng.h>

DoubleSectionChooser::DoubleSectionChooser(int nlevel, int npart,
        Paths &paths, Action &action, DoubleAction& doubleAction,
        const BeadFactory &beadFactory)
:   SectionChooser(nlevel, npart, paths, action, beadFactory),
    doubleAction(doubleAction),
    beads1(SectionChooser::beads),
    beads2(beadFactory.getNewBeads(npart,(int)pow(2,nlevel)+1)),
    permutation1(SectionChooser::permutation),
    permutation2(new Permutation(npart)) {
}

DoubleSectionChooser::~DoubleSectionChooser() {
    activateSection(1);
    delete &beads2;
    delete &permutation2;
}

void DoubleSectionChooser::run() {
  //double x; gsl_qrng_get(qrng,&x);
  double x=RandomNumGenerator::getRand()*(1-1e-8);
  int ilo=paths->getLowestOwnedSlice(true)-1;
  int ihi=paths->getHighestSampledSlice(beads->getNSlice()-1,true);
  iFirstSlice1=ilo+(int)((ihi+1-ilo)*x);
  if (iFirstSlice1>ihi) iFirstSlice1=ihi;
  if (RandomNumGenerator::getRand()>=0.5) {
    iFirstSlice1=(iFirstSlice1+paths->getNSlice()/2+1)%paths->getNSlice()-1;
  }
  //iFirstSlice1=(int)(paths.getNSlice()*x);
  iFirstSlice2=(iFirstSlice1+paths->getNSlice()/2+1)%paths->getNSlice()-1;
  activateSection(1);
  // Copy coordinates from allBeads to section Beads.
  paths->getBeads(iFirstSlice1,*beads1);
  paths->getBeads(iFirstSlice2,*beads2);
  permutation1->reset(); 
  permutation2->reset(); 
  // Initialize the action.
  action->initialize(*this);
  doubleAction.initialize(*this);
  // Run the sampling algorithm.
  CompositeAlgorithm::run();

  // Copy moved coordinates from sectionBeads to allBeads.
  paths->putDoubleBeads(iFirstSlice1,*beads1,*permutation1,
                       iFirstSlice2,*beads2,*permutation2);
  // Refresh the buffer slices.
  paths->setBuffers();
}

void DoubleSectionChooser::activateSection(const int i) {
  if (i==1) {
    beads=beads1; permutation=permutation1; iFirstSlice=iFirstSlice1;
  } else {
    beads=beads2; permutation=permutation2; iFirstSlice=iFirstSlice2;
  }
}
