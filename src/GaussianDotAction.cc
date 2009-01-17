// $Id: GaussianDotAction.cc,v 1.2 2007/12/28 22:27:04 lzhang51 Exp $
/*  Copyright (C) 2007 John B. Shumway, Jr.

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
#include "GaussianDotAction.h"
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"

GaussianDotAction::GaussianDotAction(const double v0, const double alpha,
    const Vec &center, const SimulationInfo& simInfo) 
  : v0(v0), alpha(alpha), tau(simInfo.getTau()), center(center)  {
  std::cout << "Gaussian action, alpha=" << alpha
            << ", v0= " << v0 << ", center=" << center << std::endl;
}

double GaussianDotAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  const int nTot=sectionBeads.getNPart();
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      // Add action for moving beads.
      Vec delta=movingBeads(iMoving,islice)-center;
      cell.pbc(delta);
      deltaAction+=v0*exp(-alpha*dot(delta,delta))*tau*nStride;
      // Subtract action for old beads.
      delta=sectionBeads(index(iMoving),islice)-center;
      cell.pbc(delta);
      deltaAction-=v0*exp(-alpha*dot(delta,delta))*tau*nStride;
    }
  }
  return deltaAction;
}

double GaussianDotAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void GaussianDotAction::getBeadAction(const Paths& paths, int ipart, int islice,
         double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  Vec delta=paths(ipart,islice)-center;
  paths.getSuperCell().pbc(delta);
  double v=v0*exp(-alpha*dot(delta,delta));
  utau+=v;
  //fm=delta; fm*=v*denom*denom;
  //fp=delta; fp*=v*denom*denom;
  u=utau*tau;
}
