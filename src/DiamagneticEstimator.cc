//$Id: DiamagneticEstimator.cc 22 2009-03-06 13:52:07Z john.shumwayjr $
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
#include "DiamagneticEstimator.h"
#include "SimulationInfo.h"
#include "Species.h"
#include "Paths.h"
#include <blitz/tinyvec-et.h>
#include <iostream>

DiamagneticEstimator::DiamagneticEstimator(const SimulationInfo& simInfo, const double temperature)
  : ScalarEstimator("<area**2>"),
    sus(0), area(0), norm(0),temperature(temperature),q(simInfo.getNPart()){
    for (int i=0; i<q.size(); ++i) q(i)=simInfo.getPartSpecies(i).charge;
    std::cout << "Area Squared Estimator Testing...." << std::endl;
}

void DiamagneticEstimator::initCalc(const int nslice, const int firstSlice) {
  area=0;
}

void DiamagneticEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
  area += 0.5*(start[0]*(end[1]-start[1])-start[1]*(end[0]-start[0]));
}

void DiamagneticEstimator::endCalc(const int nslice) {
//calculates <area*area>, below converts to units of metres^3  
sus+=((temperature*4.3597E-18*2.5664E-38*1.26E-6)*(area*area)*7.3116E-42)/1.1122E-68; 
//See, Pollock and Runge: Study of diamagnetic response. J Chem Phys Vol 96, No 1. 1 Jan 1992
//Equation 7 (above equation) & Equation 26 (which includes time-step error)
norm+=1;
}
