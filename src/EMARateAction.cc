// $Id$
/*  Copyright (C) 2010-11 John B. Shumway, Jr.

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
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "DisplaceMoveSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"
#include "SectionChooser.h"

EMARateAction::EMARateAction(const SimulationInfo& simInfo,
  const Species& species1, const Species& species2, double C) 
  : invTau(1./simInfo.getTau()), species1(species1), species2(species2),
    index1(species1.ifirst), index2(species2.ifirst), C(C),
    nPathSlice(simInfo.getNSlice()) {
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
  // Only evaluate if we are aligned with slice 0 in the middle.
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();

  const int iFirstSlice = sampler.getSectionChooser().getFirstSliceIndex();
  if (iFirstSlice + nSlice/2 != nPathSlice) return 0.;
  
  // For now, we'll assume that the only two radiating particles are being moved.
  int iMoving1 = index(0);
  int iMoving2 = index(1);

  double oldDiagAction = 0.;
  double oldRadAction = 0.;
  double newDiagAction = 0.;
  double newRadAction = 0.;

  const Vec inv2Sigma21 = 0.25*mass1*invTau/nStride;
  const Vec inv2Sigma22 = 0.25*mass2*invTau/nStride;
  Vec rePrev = movingBeads(0,0);
  Vec rhPrev = movingBeads(1,0);
  Vec reRadPrev = rePrev;
  Vec rhRadPrev = rhPrev;
  Vec rePrevOld = movingBeads(iMoving1,0);
  Vec rhPrevOld = movingBeads(iMoving2,0);
  Vec reRadPrevOld = rePrevOld;
  Vec rhRadPrevOld = rhPrevOld;
  for (int islice=nStride; islice<nSlice; islice+=nStride) {

    // Calculate action for moving beads.
    Vec re = movingBeads(0,islice);
    Vec rh = movingBeads(1,islice);

    Vec reRad = re;
    Vec rhRad = (islice==nSlice/2)?re:rh;

    Vec delta=re-rePrev; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      newDiagAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
    }
    delta=rh-rhPrev; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      newDiagAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
    }
    delta=reRad-reRadPrev; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      newRadAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
    }
    delta=rhRad-rhRadPrev; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      newRadAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
    }
    // Calculate action for old beads.
    Vec reOld = sectionBeads(iMoving1,islice);
    Vec rhOld = sectionBeads(iMoving2,islice);

    Vec reRadOld = re;
    Vec rhRadOld = (islice==nSlice/2)?re:rh;

    delta=reOld-rePrevOld; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      oldDiagAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
    }
    delta=rhOld-rhPrevOld; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      oldDiagAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
    }
    delta=reRadOld-reRadPrevOld; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      oldRadAction+=delta[idim]*delta[idim]*inv2Sigma21[idim];
    }
    delta=rhRadOld-rhRadPrevOld; cell.pbc(delta);
    for (int idim=0;idim<NDIM;++idim) {
      oldRadAction+=delta[idim]*delta[idim]*inv2Sigma22[idim];
    }

    // Set the previous positions.
    rePrev = re;
    rhPrev = rh;
    reRadPrev = (islice==nSlice/2)?rhPrev:rePrev;
    rhRadPrev = rhPrev;
    rePrevOld = reOld;
    rhPrevOld = rhOld;
    reRadPrevOld = (islice==nSlice/2)?rhPrevOld:rePrevOld;
    rhRadPrevOld = rhPrevOld;
  } 

  //double oldAction = -log(1+C*exp(oldRadAction-oldDiagAction));
  //double newAction = -log(1+C*exp(newRadAction-newDiagAction));

  double oldAction = -log(1+C*exp(-oldRadAction+oldDiagAction));
  double newAction = -log(1+C*exp(-newRadAction+newDiagAction));
  if (C > 0.) {
    if (log(C)-oldRadAction+oldDiagAction > 40) {
      oldAction = -log(C)+oldRadAction-oldDiagAction;
    }
    if (log(C)-newRadAction+newDiagAction > 40) {
      newAction = -log(C)+newRadAction-newDiagAction;
    }
  }

  double deltaAction = newAction-oldAction;

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
