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
#include "SHOAction.h"
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "Species.h"

SHOAction::SHOAction(const double tau, const double omega, const double mass,
                     const int ndim, const Species &species) 
  : tau(tau), omega(omega), mass(mass), ndim(ndim),
    ifirst(species.ifirst), npart(species.count) {
  std::cout << "SHOAction with omega=" << omega << std::endl;
}

double SHOAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  double wt=omega*tau*nStride;
  double sinhwt=sinh(wt);
  double coshwt=cosh(wt);
  double m=mass;
  double div1=1.0/(2.0*sinhwt);
  double div2=1.0/(2.0*tau*nStride);
  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      if (i<ifirst || i>ifirst+npart) continue;
      // Add action for moving beads.
      Vec r1=movingBeads(iMoving,islice); 
      Vec r0=movingBeads(iMoving,islice-nStride); 
      Vec delta=r1-r0;
      double r0r0=0, r0r1=0, r1r1=0, delta2=0;
      for (int idim=(NDIM-ndim); idim<NDIM; ++idim) {
        r0r0+=r0[idim]*r0[idim];
        r0r1+=r0[idim]*r1[idim];
        r1r1+=r1[idim]*r1[idim];
        delta2+=delta[idim]*delta[idim];
      }
      deltaAction+=m*omega*((r1r1+r0r0)*coshwt-2.0*r0r1)*div1-m*delta2*div2;
      // Subtract action for old beads.
      r1=sectionBeads(i,islice); r0=sectionBeads(i,islice-nStride);
      delta=r1-r0;
      r0r0=0; r0r1=0; r1r1=0; delta2=0;
      for (int idim=(NDIM-ndim); idim<NDIM; ++idim) {
        r0r0+=r0[idim]*r0[idim];
        r0r1+=r0[idim]*r1[idim];
        r1r1+=r1[idim]*r1[idim];
        delta2+=delta[idim]*delta[idim];
      }
      deltaAction-=m*omega*((r1r1+r0r0)*coshwt-2.0*r0r1)*div1-m*delta2*div2;
    }
  }
  return deltaAction;
}

double SHOAction::getTotalAction(const Paths& paths, const int level) const {
  return 0;
}

void SHOAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  u=utau=0; fm=0; fp=0;
  if (ipart<ifirst || ipart>ifirst+npart) return;
  double wt=omega*tau;
  double sinhwt=sinh(wt);
  double coshwt=cosh(wt);
  double cschwt=1.0/sinhwt;
  double cothwt=coshwt*cschwt;
  double tin=1.0/tau;
  double m=mass;
  Vec r1=paths(ipart,islice);
  Vec r0=paths(ipart,islice,-1);
  Vec r2=paths(ipart,islice,1);
  Vec delta = paths.delta(ipart,islice,-1);
  double r0r0=0, r0r1=0, r1r1=0, delta2=0;
  for (int idim=(NDIM-ndim); idim<NDIM; ++idim) {
    r0r0+=r0[idim]*r0[idim];
    r0r1+=r0[idim]*r1[idim];
    r1r1+=r1[idim]*r1[idim];
    delta2+=delta[idim]*delta[idim];
  }
  fm-= (m*omega*(r1*coshwt-r0)*cschwt-m*delta*tin);
  utau = 0.5*ndim*omega*cothwt
        +0.5*m*omega*omega
          *( (r1r1+r0r0)*(1-cothwt*cothwt)
             +2.0*r0r1*coshwt*cschwt*cschwt)
        +0.5*m*delta2*tin*tin-0.5*ndim*tin;
  u=0.5*ndim*log(2*3.14159265358979*sinhwt/(m*omega))+m*omega*((r1r1+r0r0)
	  *coshwt-2.0*r0r1) -0.5*log(2.0*3.14159265358979*tau/m)
	  -0.5*m*delta2/tau;
  delta=paths.delta(ipart,islice,1);
  fp-= (m*omega*(r1*coshwt-r2)*cschwt-m*delta*tin);
  for (int idim=0; idim<(NDIM-ndim); ++idim) {fm[idim]=0; fp[idim]=0;}
}
