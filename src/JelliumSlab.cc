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
#include "JelliumSlab.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include <fstream>

JelliumSlab::JelliumSlab(const SimulationInfo& simInfo, const double qtot,
  const double slabWidth, const bool isSlabInMiddle) 
  : tau(simInfo.getTau()), q(simInfo.getSpecies(0).charge),
    lx(simInfo.getSuperCell()->a[0])  {
  double crossSection=product(simInfo.getSuperCell()->a)/lx;
  double rhoSlab=qtot/crossSection/slabWidth;
  double rhoBackground=-qtot/crossSection/lx;
  if (isSlabInMiddle) {
    w=slabWidth;
    rhoMid=rhoSlab+rhoBackground;
    rhoEdge=rhoBackground;
  } else {
    w=lx-slabWidth;
    rhoMid=rhoBackground;
    rhoEdge=rhoSlab+rhoBackground;
  }
  double phiSum=4*PI/(3.*lx)*(rhoMid*pow(0.5*w,3)+rhoEdge*pow(0.5*(lx-w),3));
  double phiDiff=2*PI*(rhoMid*pow(0.5*w,2)-rhoEdge*pow(0.5*(lx-w),2));
  phiMid=(phiSum/(lx-w)+phiDiff)/(1+w/(lx-w));
  phiEdge=(phiSum/w-phiDiff)/(1+(lx-w)/w);

  std::ofstream file("slab.dat");
  for (int i=-50; i<50; ++i) {
    double x=lx*i*0.01;
    file << x << " " << phi(x) << std::endl;
  }
}

double JelliumSlab::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  double ktstride = 0.5*tau*nStride;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads.
      Vec delta=movingBeads(iMoving,islice);
      cell.pbc(delta);
      deltaAction+=q*phi(delta[0])*ktstride;
      // Subtract action for old beads.
      delta=sectionBeads(i,islice);
      cell.pbc(delta);
      deltaAction-=q*phi(delta[0])*ktstride;
    }
  }
  return deltaAction;
}

double JelliumSlab::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void JelliumSlab::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  Vec delta=paths(ipart,islice);
  utau=q*phi(delta[0]);
  u=utau*tau;
  fm=-0.5*tau*delta;
  fp=-0.5*tau*delta;
}

const double JelliumSlab::PI=acos(-1.0);
