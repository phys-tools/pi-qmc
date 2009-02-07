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
#include "VVCorrelationEstimator.h"

VVCorrelationEstimator::VVCorrelationEstimator(
  const SimulationInfo& simInfo)
  : ScalarEstimator("VVCorrelation"),
    v(simInfo.getNSlice()), vdot(0), vtot(0), vnorm(0), 
    cell(*simInfo.getSuperCell()) {
}

void VVCorrelationEstimator::initCalc(const int nslice, const int firstSlice) {
  v=0.;
}

void VVCorrelationEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  Vec diff=end-start;
  cell.pbc(diff);
  v(islice)+=diff;
}


void VVCorrelationEstimator::endCalc(const int nslice) {
  vdot=0.;
  int nsliceon2=nslice/2; //what about odd number of slices.
  for (int i=0; i<nsliceon2; i++)
    vdot+=blitz::dot(v(i),v(i+nsliceon2));
  vdot/=nsliceon2;
  vtot+=vdot;
  vnorm+=1;
}
