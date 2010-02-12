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
#include "GateAction.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "SuperCell.h"
#include "Paths.h"
#include "SimulationInfo.h"
#include "Species.h"
#include <stdio.h>
#include <math.h>

/*
This is a potential mimicking the squre gate potential in experiments.
The potential is GVolt * ((tanh(sx(x + xwidth/2)) - tanh(sx(x - xwidth/2))))
                       * ((tanh(sy(y + ywidth/2)) - tanh(sy(y - ywidth/2))))
GVolt is the potential of the top gate in energy unit.
s is a parameter controling the flatness of the top of the potential.
offset is a parameter representing the offset w.r.t. x-axis and y-axis.

The action is not normalized to GVolt !!
Virial -- not works

*/

GateAction::GateAction(const SimulationInfo &simInfo, const double GVolt,
 const double sx, const double sy, const double xwidth, const double ywidth, 
 const double xoffset, const double yoffset, const Species &species)
  : tau(simInfo.getTau()), GVolt(GVolt), sx(sx), sy(sy),  xwidth(xwidth), 
    ywidth(ywidth), xoffset(xoffset), yoffset(yoffset), ifirst(species.ifirst),
    npart(species.count) {
}

double GateAction::getActionDifference(const MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2, level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex();
  const int nMoving=index.size();
  double deltaAction=0;
  double ktstride = tau*nStride;
  double xspread = xwidth / 2;
  double yspread = ywidth / 2;
  double normalConst = (tanh(sx * xspread) - tanh(sx * (-xspread)))
                     * (tanh(sy * yspread) - tanh(sy * (-yspread)));
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      if (i<ifirst || i>+ifirst+npart) continue;
      // Add action for moving beads.
      Vec delta=movingBeads(iMoving, islice);
      cell.pbc(delta);
      double x=0, y=0;
      x = delta[0];
      y = delta[1];
      deltaAction += GVolt * (tanh(sx * (x - xoffset + xspread)) - 
                              tanh(sx * (x - xoffset - xspread))) *
                             (tanh(sy * (y - yoffset + yspread)) - 
                              tanh(sy * (y - yoffset - yspread))) * ktstride 
                             / normalConst;
      // Subtract action for old beads.
      delta = sectionBeads(i, islice);
      cell.pbc(delta);
      x = delta[0];
      y = delta[1];
      deltaAction -= GVolt * (tanh(sx * (x - xoffset + xspread)) - 
                              tanh(sx * (x - xoffset - xspread))) *
                             (tanh(sy * (y - yoffset + yspread)) - 
                              tanh(sy * (y - yoffset - yspread))) * ktstride 
                             / normalConst;
    }
  }
  return deltaAction;
}

double GateAction::getTotalAction(const Paths& paths, const int level) const {
  return 0;
}

void GateAction::getBeadAction(const Paths& paths, int ipart, int islice, double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
  Vec delta = paths(ipart,islice);
  fm = 0; fp = 0; ulambda = 0;
  if (ipart < ifirst || ipart >= ifirst + npart) return;
  double x = delta[0];
  double y = delta[1];
  double xspread = xwidth / 2;
  double yspread = ywidth / 2;
  double normalConst = (tanh(sx * xspread) - tanh(sx * (-xspread)))
                     * (tanh(sy * yspread) - tanh(sy * (-yspread)));
  utau = GVolt * (tanh(sx * (x - xoffset + xspread)) - 
                          tanh(sx * (x - xoffset - xspread))) *
                         (tanh(sy * (y - yoffset + yspread)) - 
                          tanh(sy * (y - yoffset - yspread))) ;
  u = utau * tau / normalConst;
}
