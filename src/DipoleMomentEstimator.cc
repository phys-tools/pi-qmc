//$Id$
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
#include "DipoleMomentEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

DipoleMomentEstimator::DipoleMomentEstimator(const SimulationInfo& simInfo, const int index)
  : ScalarEstimator(index==0 ? "dipole_moment_x" 
                  : index==1 ? "dipole_moment_y"
                             : "dipole_moment_z", 
                    "scalar-length/dipole","",1.,0.),
    dipole(0), dtot(0), dnorm(0), q(simInfo.getNPart()), index(index) {
    for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;
    std::cout << "Dipole Moment Estimator component = "
              << (index==0 ? "x" : index==1 ? "y" : "z") << std::endl;
}

void DipoleMomentEstimator::initCalc(const int nslice, const int firstSlice) {
  dipole=0;
}

void DipoleMomentEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  dipole+=q(ipart)*end(index);
}

void DipoleMomentEstimator::endCalc(const int nslice) {
  dipole/=nslice;
  dtot+=dipole; dnorm+=1;
}
