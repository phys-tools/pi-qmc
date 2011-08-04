// $Id$
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#include "OpticalLatticeAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"

OpticalLatticeAction::OpticalLatticeAction(const Vec v0, const Vec length,
  const Vec max, const SimulationInfo &simInfo)
  : tau(simInfo.getTau()), v0(v0), piOverL(3.14159265354/length), max(max) {
}

double OpticalLatticeAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  double ktstride = tau*nStride;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads.
      Vec r=movingBeads(iMoving,islice);
      cell.pbc(r);
      double pot=0;
      for (int i=0; i<NDIM; ++i) {
        if (fabs(r[i])<max[i]) {
          double sinkx = sin(piOverL[i]*r[i]);
          pot += v0[i]*(sinkx*sinkx-0.5);
        } else {
          pot += 0.5*fabs(v0[i]);
        }
      }
      deltaAction += pot*ktstride;
      // Subtract action for old beads.
      r=sectionBeads(i,islice);
      cell.pbc(r);
      pot=0;
      for (int i=0; i<NDIM; ++i) {
        if (fabs(r[i])<max[i]) {
          double sinkx = sin(piOverL[i]*r[i]);
          pot += v0[i]*(sinkx*sinkx-0.5);
        } else {
          pot += 0.5*fabs(v0[i]);
        }
      }
      deltaAction -= pot*ktstride;
    }
  }
  return deltaAction;
}

double OpticalLatticeAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void OpticalLatticeAction::getBeadAction(const Paths& paths, int ipart, 
    int islice, double& u, double& utau, double& ulambda, 
    Vec &fm, Vec &fp) const {
  Vec r=paths(ipart,islice);
  for (int i=0; i<NDIM; ++i) {
    if (fabs(r[i])<max[i]) {
      double sinkx = sin(piOverL[i]*r[i]);
      utau += v0[i]*(sinkx*sinkx-0.5);
    } else {
      utau += 0.5*fabs(v0[i]);
    }
  }
  u=utau*tau;
}
