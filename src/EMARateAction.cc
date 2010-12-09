// $Id$
/*  Copyright (C) 2010 John B. Shumway, Jr.

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
#include "EMARateAction.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "DisplaceMoveSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"

EMARateAction::EMARateAction(const SimulationInfo& simInfo,
  const Species& species1, const Species& species2, double C) 
  : invTau(1./simInfo.getTau()), species1(species1), species2(species2),
    index1(species1.ifirst), index2(species2.ifirst), C(C)  {
  if (species1.anMass) {
    mass1 = *species1.anMass;
  } else {
    mass1 = species1.mass;
  }
  if (species2.anMass) {
    mass2 = *species2.anMass;
  } else {
    mass2 = species2.mass;
  }
}

EMARateAction::~EMARateAction() {
}

double EMARateAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  //const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  //const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  //const SuperCell& cell=sampler.getSuperCell();
  //const int nStride=(int)pow(2,level);
  //const int nSlice=sectionBeads.getNSlice();
  //const IArray& index=sampler.getMovingIndex(); 
  //const int nMoving=index.size();
  double deltaAction=0;
  if (level==0) {
  }
  return deltaAction;
}

double EMARateAction::getActionDifference(const Paths &paths, 
   const VArray &displacement, int nmoving, const IArray &movingIndex, 
   int iFirstSlice, int nslice) {
 return 0; //No change in action for uniform displacements of particles.
}

double EMARateAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void EMARateAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double &u, double &utau, double &ulambda, Vec &fm, Vec &fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
}
