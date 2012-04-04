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
#include "GaussianAction.h"
#include "sampler/SectionSamplerInterface.h"
#include "Beads.h"
#include "Paths.h"
#include "util/SuperCell.h"
#include "SimulationInfo.h"

GaussianAction::GaussianAction(const double v0, const double alpha,
                             const SimulationInfo& simInfo) 
  : v0(v0), alpha(alpha), tau(simInfo.getTau())  {
  std::cout << "Gaussian action, alpha=" << alpha
            << ", v0= " << v0 << std::endl;
}

double GaussianAction::getActionDifference(const SectionSamplerInterface& sampler,
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
      const int i=index(iMoving);
      for (int j=0; j<nTot; ++j) {
        if (i==j) continue;
        // Add action for moving beads.
        Vec delta=movingBeads(iMoving,islice);
        delta-=sectionBeads(j,islice);
        cell.pbc(delta);
        deltaAction+=v0*exp(-alpha*sqrt(dot(delta,delta)))*tau*nStride;
        // Subtract action for old beads.
        delta=sectionBeads(i,islice);
        delta-=sectionBeads(j,islice);
        cell.pbc(delta);
        deltaAction-=v0*exp(-alpha*sqrt(dot(delta,delta)))*tau*nStride;
      }
    }
  }
  return deltaAction;
}

double GaussianAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void GaussianAction::getBeadAction(const Paths& paths, int ipart, int islice,
         double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  for (int j=0; j<paths.getNPart(); ++j) {
    if (ipart==j) continue;
    Vec delta=paths(ipart,islice);
    delta-=paths(j,islice);
    paths.getSuperCell().pbc(delta);
    double v=v0*exp(-alpha*sqrt(dot(delta,delta)));
    utau+=v;
    //fm=delta; fm*=v*denom*denom;
    //fp=delta; fp*=v*denom*denom;
  }
  u=utau*tau;
}
