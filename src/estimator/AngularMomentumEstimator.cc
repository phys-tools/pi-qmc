#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "AngularMomentumEstimator.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "fixednode/PhaseModel.h"

AngularMomentumEstimator::AngularMomentumEstimator(
  const SimulationInfo& simInfo, PhaseModel* phase)
  : ScalarEstimator("angular_momentum","scalar/angular-momentum","",1.,0.),
    angM(0), angMtot(0), angMnorm(0),
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
