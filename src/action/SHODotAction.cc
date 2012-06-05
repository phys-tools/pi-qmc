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
#include "SHODotAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "sampler/SectionSamplerInterface.h"
#include "Beads.h"
#include "util/SuperCell.h"
#include "Paths.h"
#include "Species.h"
#include "SimulationInfo.h"
#include "util/RandomNumGenerator.h"

SHODotAction::SHODotAction(const SimulationInfo &simInfo, double t, double v0, 
double omega, double z, const Species &species, 
  const bool semiclassical, const int numds)
  : tau(simInfo.getTau()), t(t), v0(v0), k(species.mass*omega*omega), z(z),
    ifirst(species.ifirst), npart(species.count), semiclassical(semiclassical), numds(numds), cell(*simInfo.getSuperCell()), nslice(simInfo.getNSlice()), 
    nds1(1.0/numds) {
  std::cout << "SHODotAction: " << species << std::endl;
  std::cout << "thick,v0,mw2,z: " << t << "," << v0 << "," << k 
  << ", " << z << std::endl;
  if(semiclassical) {
    std::cout <<"Using Semiclassical improvement for SHODotAction with numds "<<numds<<std::endl;
  }
}

double SHODotAction::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
#if NDIM==3
  double ktstride = 0.5*k*tau*nStride;
#endif
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    int islicePrime = islice-nStride;
    if( islicePrime<0 )
      islicePrime+=nslice;
    else if( islicePrime>=nslice )
      islicePrime-=nslice;

    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      if (i<ifirst || i>=ifirst+npart) continue;
      // Add action for moving beads.
      Vec position = movingBeads(iMoving,islice);
      cell.pbc(position);
#if NDIM==3
      deltaAction+=ktstride*(position[0]*position[0]+position[1]*position[1]);
      if(semiclassical)
        deltaAction+=tau*getSemiClassicalAction(movingBeads(iMoving,islicePrime)[2],position[2], nStride);
      else
        deltaAction+=(fabs(position[2]-z)>0.5*t)?v0*tau*nStride:0;
#else
      deltaAction+=(fabs(position[0]-z)>0.5*t)?v0*tau*nStride:0;
#endif
      // Subtract action for old beads.
      position=sectionBeads(i,islice);
      cell.pbc(position);
#if NDIM==3
      deltaAction-=ktstride*(position[0]*position[0]+position[1]*position[1]);
      if(semiclassical)
        deltaAction-=tau*getSemiClassicalAction(movingBeads(iMoving,islicePrime)[2], position[2], nStride);
      else
        deltaAction-=(fabs(position[2]-z)>0.5*t)?v0*tau*nStride:0;
#else
      deltaAction-=(fabs(position[0]-z)>0.5*t)?v0*tau*nStride:0;
#endif
    }
  }
  return deltaAction;
}

double SHODotAction::getActionDifference(const Paths &paths, const VArray &displacement, int nmoving, const IArray &movingIndex, int iFirstSlice, int iLastSlice) {
  double deltaAction = 0;
#if NDIM==3
  double kt = 0.5*k*tau;
#endif
  const SuperCell& cell = paths.getSuperCell();
  for (int i=0; i<nmoving; ++i) {
    int ipart = movingIndex(i);
    if (ipart<ifirst || ipart>=ifirst+npart) break;
    for (int islice=iFirstSlice; islice<=iLastSlice; ++islice) {
      Vec position = paths(ipart,islice);
      cell.pbc(position);

#if NDIM==3
      deltaAction-=kt*(position[0]*position[0]+position[1]*position[1]);
      deltaAction-=(fabs(position[2]-z)>0.5*t)?v0*tau:0;
#else
      deltaAction-=(fabs(position[0]-z)>0.5*t)?v0*tau:0;
#endif

      position += displacement(i);
      cell.pbc(position);

#if NDIM==3
      deltaAction+=kt*(position[0]*position[0]+position[1]*position[1]);
      deltaAction+=(fabs(position[2]-z)>0.5*t)?v0*tau:0;
#else
      deltaAction+=(fabs(position[0]-z)>0.5*t)?v0*tau:0;
#endif
    }
  }
  return deltaAction;
}

double SHODotAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void SHODotAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  double z1,z2;
  if (ipart<ifirst || ipart>=ifirst+npart) return;
  Vec position=paths(ipart,islice);
  if(semiclassical) {
    z1 = position[2];
    z2 = paths(ipart,islice,-1)[2];
  }
#if NDIM==3
//  position[2] -= z;
//  cell.pbc(position);
  utau=0.5*k*(position[0]*position[0]+position[1]*position[1]);
  if(semiclassical)
    utau+=getSemiClassicalAction(z1, z2, 1);
  else
    utau+=(fabs(position[2]-z)>0.5*t)?v0:0;
#else
  utau=(fabs(position[0]-z)>0.5*t)?v0:0;
#endif
  u=utau*tau;
}

//Semiclassical improvement to the action based on (4.13) of Ceperley RMP
//we leave out a factor of tau to be multiplied later
//FIXME: Probably breaks workers
double SHODotAction::getSemiClassicalAction(double zinitial, double zfinal,
  int nStride) const {

  double deltaAction=0.0; 
  double deltaz, zcurrent;
#if NDIM==3
  deltaz = cell.pbc(zfinal-zinitial,2);
  deltaz*=nds1;
  for (int ids=0;ids<=numds;ids++) {
    zcurrent=cell.pbc(zinitial+(ids*deltaz),2);
    deltaAction+=(fabs(zcurrent-z)>0.5*t)?1:0;
  }
  deltaAction*=(v0*nStride*nds1);
#endif
  return deltaAction;
}
