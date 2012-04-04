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
#include "sampler/SectionSamplerInterface.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include "SuperCell.h"
#include <blitz/tinyvec-et.h>
#include <iostream>
#include <fstream>

EFieldAction::EFieldAction(const SimulationInfo& simInfo,
    double strength, double center, double width, int idir)
  : tau(simInfo.getTau()), q(simInfo.getNPart()), strength(strength), 
    center(0.), halfwidth(0.5*width), idir(idir),
    cell(*simInfo.getSuperCell()),
    slopeIn(-strength), slopeOut(strength*width/(cell.a[idir]-width)),
    intercept(0.5*cell.a[idir]*slopeOut) {
  std::cout << (idir==0 ? "EFieldAction x" 
              : idir==1 ? "EFieldAction y"
                        : "EFieldAction z") << std::endl;
  std::cout << "strength = " << strength << ", ";
  std::cout << "center = " << center << ", ";
  std::cout << "width = " << width << ", ";
  std::cout << (idir==0 ? "X" : idir==1 ? "Y" : "Z") << std::endl;
  std::cout << slopeIn << ", " << slopeOut << ", " << intercept << std::endl;
  //std::cout << " = " << cell.a[idir] << std::endl;
  for(int i=0; i<q.size(); i++) {
    q(i)=simInfo.getPartSpecies(i).charge;
  }
  std::ofstream file("efield.out");
  for (int i=0; i<=1000; ++i) {
    double x=0.001*(i-500)*cell.a[idir];
    file << x << " " << v(x) << std::endl;
  }
}

double EFieldAction::getActionDifference(const SectionSamplerInterface& sampler,
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
      Vec r2=movingBeads(iMoving, islice)-center;
      cell.pbc(r2);
      // old beads.
      Vec r1=sectionBeads(i, islice)-center;
      cell.pbc(r1);
      deltaAction+=tau*nStride*q(i)*(v(r2(idir))-v(r1(idir)));
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
    Vec delta=paths(ipart,islice)-center;
    cell.pbc(delta);
    utau=q(ipart)*v(delta(idir));
    u=utau*tau;
}

double EFieldAction::v(const double z) const {
  // Calculate all three cases to avoid an if statement.
  double vmid = slopeIn*z;
  double vhi =  slopeOut*z-intercept;
  double vlo =  slopeOut*z+intercept;
  const double lookup_table[] = {vmid,vhi,vlo};
  return lookup_table[(z>halfwidth) + (z<-halfwidth)*2];
}

double EFieldAction::getActionDifference(const Paths &paths,
       const VArray &displacement, int nmoving, const IArray &movingIndex,
       int iFirstSlice, int iLastSlice) {

  const SuperCell& cell=paths.getSuperCell();
  double deltaAction=0;
  for (int iMoving=0; iMoving<nmoving; ++iMoving) {
    int ipart = movingIndex(iMoving);

    for (int islice=iFirstSlice; islice<=iLastSlice; ++islice) {
      // Subtract action for old positions. (Evaluate v at midpoint)
      
      Vec r = paths(ipart,islice,-1)-center;
      cell.pbc(r);
      Vec delta = paths(ipart,islice);
      delta-=r; cell.pbc(delta);delta*=0.5; r+=delta;
      cell.pbc(r);
      deltaAction -= v(r(idir))*q(ipart)*tau;

      // Add action for displaced position.
      r += displacement(iMoving);
      cell.pbc(r);
      deltaAction += v(r(idir))*q(ipart)*tau;

    }
  }
  return deltaAction;
}
