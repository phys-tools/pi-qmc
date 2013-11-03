#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GroundStateWFNodes.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/Species.h"
#include "util/SuperCell.h"

GroundStateWFNodes::GroundStateWFNodes(const Species& s, const Species& rs)
  : npart(s.count),ifirst(s.ifirst), iref(rs.ifirst) {
}

double GroundStateWFNodes::getActionDifference(const SectionSamplerInterface& sampler,
                                               int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride = 1 << level;
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  // Figure out which particles are moving.
  isMoving0=false; isMoving1=false; isRefMoving=false;
  for (int iMoving=0; iMoving<nMoving; ++iMoving) {
    const int i=index(iMoving);
    if (i==iref) {isRefMoving=true; indref=iMoving;}
    if (i==ifirst) {isMoving0=true; ind0=iMoving;}
    if (i==ifirst+1) {isMoving1=true; ind1=iMoving;}
  }
  // Calculate old determinant sign.
  Vec v0=sectionBeads(ifirst,0);
  v0-=sectionBeads(iref,0);
  cell.pbc(v0);
  Vec v1=sectionBeads(ifirst+1,0);
  v1-=sectionBeads(iref,0);
  cell.pbc(v1);
  double sign=dot(v0,v0)-dot(v1,v1);
  for (int islice=nStride; islice<nSlice; islice+=2*nStride) {
    // Now check new determinant.
    v0=(isMoving0)?movingBeads(ind0,islice):sectionBeads(ifirst,islice);
    v0-=(isRefMoving)?movingBeads(indref,islice):sectionBeads(iref,islice);
    cell.pbc(v0);
    v1=(isMoving1)?movingBeads(ind1,islice):sectionBeads(ifirst+1,islice);
    v1-=(isRefMoving)?movingBeads(indref,islice):sectionBeads(iref,islice);
    cell.pbc(v1);
    if ((dot(v0,v0)-dot(v1,v1))*sign<0) return 1.e200;
  } 
  return 0;
}
