#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "SphereAction.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/Species.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>

SphereAction::SphereAction(const double tau, const double radius,
                           const Species& species)
 : tau(tau), radius(radius), ncenter(species.count), ipart(ncenter),
   isMoving(ncenter), icm(ncenter) {
 for (int i=0; i<ncenter; ++i) ipart(i)=species.ifirst+i;
 std::cout << "SphereAction with radius=" << radius 
      << " and center(s) " << ipart << std::endl;
}

double SphereAction::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride = 1 << level;
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  isMoving=false; icm=0;
  for (int iMoving=0; iMoving<nMoving; ++iMoving) {
    for (int icenter=0; icenter<ncenter; ++icenter) {
      if (index(iMoving)==ipart(icenter)) {
        isMoving(icenter)=true; icm(icenter)=iMoving;
      }
    }
  }
  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      // Add action for moving beads.
      // Loop over centers and make sure at least one is close.
      for (int icenter=0; icenter<ncenter; ++icenter) {
        Vec center=(isMoving(icenter))
                   ?movingBeads(icm(icenter),islice)
                   :sectionBeads(ipart(icenter),islice);
        Vec delta=movingBeads(iMoving,islice)-center;
        cell.pbc(delta);
        if (dot(delta,delta)<radius*radius) break; //Close to this center.
        if (icenter==ncenter-1) return 1e200; //Outside of all spheres.
      }
    }
  }
  return 0;
}
