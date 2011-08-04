// $Id$
/*  Copyright (C) 2004-2006,2009 John B. Shumway, Jr.

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
#include "PrimSHOAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include "Species.h"

PrimSHOAction::PrimSHOAction(const double a, const double b, 
  const SimulationInfo &simInfo, int ndim, const Species &species) 
  : tau(simInfo.getTau()), a(a), b(b), ndim(ndim),
  ifirst(species.ifirst), npart(species.count){
}

double PrimSHOAction::getActionDifference(const MultiLevelSampler& sampler,
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
if (i<ifirst || i>=ifirst+npart) continue;
      // Add action for moving beads.
      Vec delta=movingBeads(iMoving,islice);
      cell.pbc(delta);
      //double x2 = dot(delta,delta);
      double x2=0;
      for (int i=NDIM-ndim; i<NDIM; ++i) x2+=delta[i]*delta[i];
      deltaAction+=(a*x2+b*x2*x2)*ktstride;
      // Subtract action for old beads.
      delta=sectionBeads(i,islice);
      cell.pbc(delta);
      //x2 = dot(delta,delta);
      x2=0;
      for (int i=NDIM-ndim; i<NDIM; ++i) x2+=delta[i]*delta[i];
      deltaAction-=(a*x2+b*x2*x2)*ktstride;
    }
  }
  return deltaAction;
}

double PrimSHOAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void PrimSHOAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
u=utau=0; fm=0; fp=0;
if (ipart<ifirst || ipart>=ifirst+npart) return;
  Vec delta=paths(ipart,islice);
  //double x2=dot(delta,delta);
  double x2=0;
  for (int i=NDIM-ndim; i<NDIM; ++i) x2+=delta[i]*delta[i];
  utau = a*x2 + b*x2*x2;
  u=utau*tau;
  fm=-0.5*(2*a+4*b*x2)*tau*delta;
  fp=-0.5*(2*a+4*b*x2)*tau*delta;
}
