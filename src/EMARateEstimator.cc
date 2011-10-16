//$Id$
/*  Copyright (C) 2011 John B. Shumway, Jr.

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
#include "EMARateEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include "SuperCell.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

EMARateEstimator::EMARateEstimator(const SimulationInfo& simInfo, double C)
  : ScalarEstimator("ema_rate","scalar/ema-rate","",1.,0.),
    dtau(simInfo.getTau()), masse(1.), massh(1.), C(C), 
    cell(*simInfo.getSuperCell()), sum(0), norm(0) {
}

void EMARateEstimator::initCalc(const int nslice, const int firstSlice) {
  actionDifference = 0.;
}

void EMARateEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  if (islice==0) {
    if (ipart==0) {
      Vec delta = start-paths(1,islice,-1);
      cell.pbc(delta);
      actionDifference += 0.5*massh*dot(delta,delta)/dtau;
    } else if (ipart==1) {
      Vec delta = start-end;
      cell.pbc(delta);
      actionDifference -= 0.5*massh*dot(delta,delta)/dtau;
    }
  } else if (islice==1) {
    if (ipart==0) {
      Vec delta = start-paths(1,islice,-1);
      cell.pbc(delta);
      actionDifference += 0.5*masse*dot(delta,delta)/dtau;
      delta = start-end;
      cell.pbc(delta);
      actionDifference -= 0.5*masse*dot(delta,delta)/dtau;
    }
  }
}

void EMARateEstimator::endCalc(const int nslice) {
//std::cout << 1./(1.+exp(-actionDifference)) << std::endl;
  if (C>0. && actionDifference-log(C) > -40) {
    sum += 1./(1.+C*exp(-actionDifference));
  } else {
    sum += 0.;
  }
  norm += 1;
}
