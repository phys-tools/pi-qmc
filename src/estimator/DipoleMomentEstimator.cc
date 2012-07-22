#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "DipoleMomentEstimator.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
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
