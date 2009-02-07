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
#include "SpringTensorAction.h"
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "Species.h"

SpringTensorAction::SpringTensorAction(const SimulationInfo& simInfo)
  : lambda(simInfo.getNPart()), tau(simInfo.getTau()) {
  for (int i=0; i<simInfo.getNPart(); ++i) {
    lambda(i)=0.5;
    lambda(i)/=(*simInfo.getPartSpecies(i).anMass);
  }
}

double SpringTensorAction::getActionDifference(
    const MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      double inv2Sigma2x = 0.25/(lambda(i)[0]*tau*nStride);
      double inv2Sigma2y = 0.25/(lambda(i)[1]*tau*nStride);
      double inv2Sigma2z = 0.25/(lambda(i)[2]*tau*nStride);
      // Add action for moving beads.
      Vec delta=movingBeads(iMoving,islice);
      delta-=movingBeads(iMoving,islice-nStride);
      cell.pbc(delta);
      deltaAction+=(delta[0]*delta[0]*inv2Sigma2x+
                    delta[1]*delta[1]*inv2Sigma2y+
                    delta[2]*delta[2]*inv2Sigma2z);
      // Subtract action for old beads.
      delta=sectionBeads(i,islice);
      delta-=sectionBeads(i,islice-nStride);
      cell.pbc(delta);
      deltaAction-=(delta[0]*delta[0]*inv2Sigma2x+
                    delta[1]*delta[1]*inv2Sigma2y+
                    delta[2]*delta[2]*inv2Sigma2z);
    }
  }
  return deltaAction;
}

double SpringTensorAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void SpringTensorAction::getBeadAction(const Paths& paths, const int ipart,
    const int islice, double &u, double &utau, double &ulambda, 
    Vec &fm, Vec &fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.; 
  Vec delta = paths.delta(ipart,islice,-1);
  utau  = 0.5/tau - delta[0]*delta[0]/(4.0*lambda(ipart)[0]*tau*tau);
  utau += 0.5/tau - delta[1]*delta[1]/(4.0*lambda(ipart)[1]*tau*tau);
  utau += 0.5/tau - delta[2]*delta[2]/(4.0*lambda(ipart)[2]*tau*tau);
}
