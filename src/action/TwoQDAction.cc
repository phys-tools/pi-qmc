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
#include "TwoQDAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "sampler/SectionSamplerInterface.h"
#include "Beads.h"
#include "util/SuperCell.h"
#include "Paths.h"

TwoQDAction::TwoQDAction(const double tau, const double omega, 
  const double mass, const double d, const double alpha)
  : tau(tau), omega(omega), mass(mass), d(d), alpha(alpha) {
  std::cout<< "TwoQDAction with omega=" << omega << ", alpha=" << alpha
	   << " and interdot distance="  << d << std::endl;
}

double TwoQDAction::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  //const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  double mw2=0.5*mass*omega*omega;
  double kt=tau*nStride;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads.
      Vec r1=movingBeads(iMoving,islice); 
      double xmd2=(r1[0]-0.5*d)*(r1[0]-0.5*d)+alpha*r1[1]*r1[1];	
      double xpd2=(r1[0]+0.5*d)*(r1[0]+0.5*d)+alpha*r1[1]*r1[1];
      if (xmd2>=xpd2){ deltaAction+= kt*mw2*xpd2 ;}
        else { deltaAction+= kt*mw2*xmd2 ;}
      // Subtract action for old beads.
      r1=sectionBeads(i,islice);
      xmd2=(r1[0]-0.5*d)*(r1[0]-0.5*d)+alpha*r1[1]*r1[1];	
      xpd2=(r1[0]+0.5*d)*(r1[0]+0.5*d)+alpha*r1[1]*r1[1];
      if (xmd2>=xpd2) {deltaAction-=kt*mw2*xpd2;}
        else {deltaAction-=kt*mw2*xmd2;}
    }
  }
  return deltaAction;
}

double TwoQDAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void TwoQDAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  Vec r1=paths(ipart,islice);
  double mw2=0.5*mass*omega*omega;
  double xmd2=(r1[0]-0.5*d)*(r1[0]-0.5*d)+alpha*r1[1]*r1[1];	
  double xpd2=(r1[0]+0.5*d)*(r1[0]+0.5*d)+alpha*r1[1]*r1[1];
  if (xmd2>=xpd2) {utau=mw2*xpd2;}
    else  {utau= mw2*xmd2;}
  u=utau*tau;
  //f=-tau*delta;
}
