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
#include "EFieldAction.h"
#include "Beads.h"
#include "MultiLevelSampler.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include "SuperCell.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

EFieldAction::EFieldAction(const SimulationInfo& simInfo,
                           const double scale, const int index)
  : tau(simInfo.getTau()), q(simInfo.getNPart()), scale(scale), 
component(index),
    cell(*simInfo.getSuperCell()), a(0.25*cell.a[component]) {
  std::cout << (component==0 ? "EFieldAction x" : component==1 ? "EFieldAction y" : "EFieldAction z") << std::endl;
  std::cout << "scale = " << scale << ", ";
  std::cout << (component==0 ? "X" : component==1 ? "Y" : "Z");
  std::cout << " = " << cell.a[component] << std::endl;
  for(int i=0; i<q.size(); i++)
    q(i)=simInfo.getPartSpecies(i).charge;
}

double EFieldAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex();
  const int nMoving=index.size();
  double deltaAction=0.;
  for(int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for(int iMoving=0; iMoving<nMoving; iMoving++) {
      const int i=index(iMoving);
      // moving beads.
      Vec r2=movingBeads(iMoving, islice);
      cell.pbc(r2);
      // old beads.
      Vec r1=sectionBeads(i, islice);
      cell.pbc(r1);
      deltaAction+=tau*nStride*q(i)*(v(r2(component))-v(r1(component)));
    }
  }
  return deltaAction;
}

double EFieldAction::getTotalAction(const Paths& paths, int level) const {
  return 0.;
}

void EFieldAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp)
    const {
    u=utau=ulambda=0; fm=0; fp=0;
    Vec delta=paths(ipart,islice);
    cell.pbc(delta);
    utau=q(ipart)*v(delta(component));
    u=utau*tau;
}

double EFieldAction::v(const double z) const {
  double v=0.;
  const double L2=cell.a[component]/2.;
  if (z<-a) {
    v=-a*(L2+z)/(L2-a);
  } else if (z<a) {
    v=z;
  } else {
    v=-a*(z-L2)/(a-L2);
  }
  return -scale*v;
/*  const double Z=cell.a[component]/2.;
  if(z>=-Z && z<-1.*Z/11.)
    v=1.1*(z/Z+1.);
  else if (z>=-1.*Z/11. && z<1.*Z/11.)
    v=-11.*z/Z;
  else if (z>=1.*Z/11. && z<Z)
    v=1.1*(z/Z-1.);
  return scale*v; */
}

