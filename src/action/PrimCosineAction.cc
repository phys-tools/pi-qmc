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
#include <stdio.h>
#endif
#include "PrimCosineAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "sampler/SectionSamplerInterface.h"
#include "Beads.h"
#include "util/SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"

PrimCosineAction::PrimCosineAction(const double a, const double b, 
  const SimulationInfo &simInfo, int ndim) 
  : tau(simInfo.getTau()), a(a), b(b), ndim(ndim) {
}

/*
1d Periodic Cosine potential, modified from PrimSHOAction
Used to test. Doesn't have an analytical answer.
Checked against a Matlab diagonalization routine available from
pgm4@hw.ac.uk (will perhaps try to put it online too)
Peter G McDonald 2009, Heriot Watt University, Scotland. 
*/

double PrimCosineAction::getActionDifference(const SectionSamplerInterface& sampler,
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
      Vec delta=movingBeads(iMoving,islice);
      cell.pbc(delta);
      //double x2 = dot(delta,delta);
      double x2=0;
      for (int i=NDIM-ndim; i<NDIM; ++i) x2+=delta[i];
      deltaAction+=(a*(1-cos(2*3.14159265*(x2/b))))*ktstride;
      // Subtract action for old beads.
      delta=sectionBeads(i,islice);
      cell.pbc(delta);

x2=delta[0];      

deltaAction-=(a*(1-cos(2*3.14159265*(x2/b))))*ktstride;

    }
  }
  return deltaAction;
}

double PrimCosineAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void PrimCosineAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  Vec delta=paths(ipart,islice);
fm=0; fp=0; ulambda=0;
  double x2=0;
  x2=delta[0];
  utau =(a*(1-cos(2*3.14159265*(x2/b))));
u=utau*tau;
fm=-a/2*(2*3.14159265/b)*sin(2*3.14159265*(x2/b))*tau;
fp=-a/2*(2*3.14159265/b)*sin(2*3.14159265*(x2/b))*tau;
}
