#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "SpinEstimator.h"
#include "action/Action.h"
#include "action/DoubleAction.h"
#include "base/SimulationInfo.h"
#include <cstdlib>
#include <blitz/tinyvec.h>

SpinEstimator::SpinEstimator(
  const SimulationInfo& simInfo, const int idim, const double gc)
  : ScalarEstimator("spin"+((idim==0)?namex:(idim==1)?namey:namez)),
    idim(idim), gc(gc), invgc(pow((1.0+gc),3)), energy(0), etot(0), enorm(0),
    l2(1/(simInfo.getSpecies(0).mass*simInfo.spinOmega)) {
}

void SpinEstimator::initCalc(const int nslice, const int firstSlice) {
  energy=0;
}

void SpinEstimator::handleLink(const Vec& start, const Vec& end,
          const int ipart, const int islice, const Paths& paths) {
#if NDIM==4 
  const Vec s=end;
#else
  const blitz::TinyVector<double,4> s
            =*reinterpret_cast<const blitz::TinyVector<double,4>*>(
               paths.getAuxBead(ipart,islice,1));
#endif
	double dots=dot(s,s)/l2;
//  energy+=dots;
	double f=1.5*invgc*exp(-gc*dots)/dots;
  switch (idim) {
    case 0:
      energy+=2*f*(-s[0]*s[2]+s[1]*s[3])/l2;
      break;
    case 1:
      energy+=2*f*(s[0]*s[1]+s[2]*s[3])/l2;
      break;
    case 2:
      energy+=f*(s[0]*s[0]-s[1]*s[1]-s[2]*s[2]+s[3]*s[3])/l2;
      break;
  }
}

void SpinEstimator::endCalc(const int nslice) {
  energy/=nslice;
  etot+=energy; enorm+=1;
}

const std::string  SpinEstimator::namex("_x"); 
const std::string  SpinEstimator::namey("_y"); 
const std::string  SpinEstimator::namez("_z"); 
