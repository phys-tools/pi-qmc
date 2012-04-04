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

DiamagneticEstimator::DiamagneticEstimator(const SimulationInfo& simInfo,
  double temperature, const std::string& unitName, double scale)
  : ScalarEstimator("chiM","scalar/diamagnetism",unitName,scale,0.),
    value(0.), area(0.), term2(0.), norm(0.),
    temperature(temperature),q(simInfo.getNPart()),
    coef(simInfo.getNPart()), tauTerm(0) {
  const double tau = simInfo.getTau();
  const int nslice = simInfo.getNSlice();
  for (int i=0; i<q.size(); ++i) {
    q(i)=simInfo.getPartSpecies(i).charge;
    double mass = simInfo.getPartSpecies(i).mass;
    coef(i) = q(i)*q(i)/(12.*mass*nslice*C*C);
    tauTerm += q(i)*q(i)*tau/(12*mass*mass*C*C);
  }
  std::cout << "Area Squared Estimator Testing...." << std::endl;
}

void DiamagneticEstimator::initCalc(const int nslice, const int firstSlice) {
  area = 0.; term2 = 0.;
}

void DiamagneticEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
#if NDIM>=2
  area += 0.5*q(ipart)*(start[0]*(end[1]-start[1])-start[1]*(end[0]-start[0]));
  term2 += coef(ipart)*((start[0]-end[0])*(start[0]-end[0])
                      + (start[1]-end[1])*(start[1]-end[1]));
#endif
}

void DiamagneticEstimator::endCalc(const int nslice) {
  /// Needs code for parallel workers!!!
  value -= temperature*(area*area)/(C*C) + term2 + tauTerm;
  norm+=1.;
}

const double DiamagneticEstimator::C = 137.035999679;
