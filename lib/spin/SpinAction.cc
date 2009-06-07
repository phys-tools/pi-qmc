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
#include "SpinAction.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "MultiLevelSampler.h"

SpinAction::SpinAction(const double tau, const double mass, 
                       const double bx, const double by, const double bz,
		       const double omega, const double gc) 
  : tau(tau), mass(mass), bx(bx), by(by), bz(bz), gc(gc),
    invgc(pow((1.0+gc),4)), omega(omega), l2(1/(mass*omega)) {
  std::cout << "Spin action with b = " << bx 
            << ", " << by << ", " << bz <<" gc=" << gc<< std::endl;
}

SpinAction::~SpinAction() {
}

double SpinAction::getActionDifference(
    const MultiLevelSampler& sampler, const int level) {
  double diff=0;
  const int nStride=(int)pow(2,level);
#if NDIM==4
  const Beads<4>& sectionBeads=sampler.getSectionBeads();
  const Beads<4>& movingBeads=sampler.getMovingBeads();
#else
  //int nAuxBeads=sampler.getSectionBeads().getAuxBeadCount();
  const Beads<4>& sectionBeads
    = *dynamic_cast<const Beads<4>*>(sampler.getSectionBeads().getAuxBeads(1));
  const Beads<4>& movingBeads
    = *dynamic_cast<const Beads<4>*>(sampler.getMovingBeads().getAuxBeads(1));
#endif
  double wt=omega*tau*nStride;
  double sinhwt=sinh(wt);
  double coshwt=cosh(wt);
  double m=mass;
  double div1=1.0/(2.0*sinhwt);
//#endif
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();

  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
     //Add new action.
      SVec s=movingBeads(iMoving,islice);
      double dots=dot(s,s);
      double f=invgc*exp(-gc*dots/l2);
      diff-=nStride*tau*( bz*f*(s[0]*s[0]-s[1]*s[1]-s[2]*s[2]+s[3]*s[3]) 
          + 2*bx*f*(-s[0]*s[2]+s[1]*s[3]) + 2*by*f*(s[0]*s[1]+s[2]*s[3]) )/l2;
//#if NDIM!=4
      SVec s0=movingBeads(iMoving,islice-nStride);
      diff+=m*omega*((dot(s,s)+dot(s0,s0))*coshwt-2.0*dot(s,s0))*div1;
//#endif
    
     //Subtract old action.
      s=sectionBeads(i,islice);
      dots=dot(s,s);
      f=invgc*exp(-gc*dots/l2);
      diff+=nStride*tau*(bz*f*( s[0]*s[0]-s[1]*s[1]-s[2]*s[2]+s[3]*s[3]) 
          + 2*bx*f*(-s[0]*s[2]+s[1]*s[3]) + 2*by*f*(s[0]*s[1]+s[2]*s[3]) )/l2;
//#if NDIM!=4
      s0=sectionBeads(i,islice-nStride);
      diff-=m*omega*((dot(s,s)+dot(s0,s0))*coshwt-2.0*dot(s,s0))*div1;
//#endif
    }
  }
  return diff;
}

double SpinAction::getTotalAction(const Paths& paths, const int level) const {
  double total=0;
  return total;
}

void SpinAction::getBeadAction(const Paths& paths, int ipart, int islice,
       double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  u=utau=ulambda=0; fm=0; fp=0;
  return;
#if NDIM==4
  int totNSlice=paths.getNSlice();
  int jslice=(islice+totNSlice/2)%totNSlice;
  SVec s=paths(0,islice);
  SVec sref=paths(0,jslice);
  SVec xref=SVec(-sref[3],-sref[2],sref[1],sref[0]);
  double ssref=dot(s,sref);
  double sxref=dot(s,xref);
  SVec gradPhi=(sref*sxref+xref*ssref)/(ssref*ssref+sxref*sxref);
  utau=dot(gradPhi,gradPhi)*0.5/mass-3*omega;
  u=utau*tau;
#endif
}

void SpinAction::acceptLastMove() {}
