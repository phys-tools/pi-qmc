// $Id: AngularMomentumEstimator.cc,v 1.2 2007/01/16 20:29:17 jshumwa Exp $
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
#include "AngularMomentumEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "PhaseModel.h"
#include "Paths.h"

AngularMomentumEstimator::AngularMomentumEstimator(
  const SimulationInfo& simInfo, PhaseModel* phase)
  : ScalarEstimator("angular_momentum"), angM(0), angMtot(0), angMnorm(0),
    phaseModel(phase), npart(simInfo.getNPart()), nslice(simInfo.getNSlice()), 
    r1(npart), r2(npart), grad(npart) {
}

AngularMomentumEstimator::~AngularMomentumEstimator() {
  delete phaseModel; phaseModel=0;
}

void AngularMomentumEstimator::initCalc(const int nslice, const int firstSlice) {
  angM=0;
}

void AngularMomentumEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
   if (ipart==0) {
     // Need to evaluate the phase at each slice.
     int jslice=(islice+nslice/2)%nslice;
     for (int i=0; i<npart; ++i) r1(i)=paths(i,islice);
     for (int i=0; i<npart; ++i) r2(i)=paths(i,jslice);
     phaseModel->evaluate(r1,r2,0);
     grad = phaseModel->getGradPhi(1);
   }
   double xdy, ydx ;
   xdy=0, ydx=0;
   xdy=paths(ipart,islice)[0]*grad(ipart)[1];
   ydx=paths(ipart,islice)[1]*grad(ipart)[0];
   angM+=xdy-ydx;
}

void AngularMomentumEstimator::endCalc(const int nslice) {
  angM/=nslice;
  angMtot+=angM; angMnorm+=1;
}
