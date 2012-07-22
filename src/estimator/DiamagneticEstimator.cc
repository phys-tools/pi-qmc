#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "DiamagneticEstimator.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "base/Paths.h"
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
