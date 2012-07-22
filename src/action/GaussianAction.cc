#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GaussianAction.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "util/SuperCell.h"

GaussianAction::GaussianAction(const double v0, const double alpha,
                             const SimulationInfo& simInfo) 
  : v0(v0), alpha(alpha), tau(simInfo.getTau())  {
  std::cout << "Gaussian action, alpha=" << alpha
            << ", v0= " << v0 << std::endl;
}

double GaussianAction::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  const int nTot=sectionBeads.getNPart();
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      for (int j=0; j<nTot; ++j) {
        if (i==j) continue;
        // Add action for moving beads.
        Vec delta=movingBeads(iMoving,islice);
        delta-=sectionBeads(j,islice);
        cell.pbc(delta);
        deltaAction+=v0*exp(-alpha*sqrt(dot(delta,delta)))*tau*nStride;
        // Subtract action for old beads.
        delta=sectionBeads(i,islice);
        delta-=sectionBeads(j,islice);
        cell.pbc(delta);
        deltaAction-=v0*exp(-alpha*sqrt(dot(delta,delta)))*tau*nStride;
      }
    }
  }
  return deltaAction;
}

double GaussianAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void GaussianAction::getBeadAction(const Paths& paths, int ipart, int islice,
         double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  for (int j=0; j<paths.getNPart(); ++j) {
    if (ipart==j) continue;
    Vec delta=paths(ipart,islice);
    delta-=paths(j,islice);
    paths.getSuperCell().pbc(delta);
    double v=v0*exp(-alpha*sqrt(dot(delta,delta)));
    utau+=v;
    //fm=delta; fm*=v*denom*denom;
    //fp=delta; fp*=v*denom*denom;
  }
  u=utau*tau;
}
