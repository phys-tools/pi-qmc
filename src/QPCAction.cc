// $Id: QPCAction.cc,v 1.5 2006/10/18 17:08:19 jshumwa Exp $
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
#include "QPCAction.h"
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"

QPCAction::QPCAction(const double tau, const double length, const double v0, 
  const double omega, const double mass)
  : tau(tau), length(length), v0(v0), omega(omega), mass(mass) {
}

double QPCAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads.
      Vec r=movingBeads(iMoving,islice); 
      //cell.pbc(r);
      deltaAction += tau*nStride*v(r[0],r[1]);
      // Subtract action for old beads.
      r=sectionBeads(i,islice);
      //cell.pbc(r);
      deltaAction -= tau*nStride*v(r[0],r[1]);
    }
  }
  return deltaAction;
}

double QPCAction::getTotalAction(const Paths& paths, const int level) const {
  return 0;
}

void QPCAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  u=utau=0; fm=0.; fp=0.;
  Vec r=paths(ipart,islice);
  utau = v(r[0],r[1]);
  u = utau*tau;
}

double QPCAction::v(double x, double y) const {
  double coshinv = 1.0/cosh(x/length);
  double vx=v0*coshinv*coshinv;
  return 0.5*vx + 0.5*mass*(omega+vx)*(omega+vx)*y*y;
}
