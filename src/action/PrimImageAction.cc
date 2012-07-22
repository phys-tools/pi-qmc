#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "PrimImageAction.h"
#include "advancer/SectionSamplerInterface.h"
#include "base/Beads.h"
#include "base/Paths.h"
#include "base/SimulationInfo.h"
#include "base/Species.h"
#include "util/SuperCell.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>

/*
This is a image charge potential, 
 
d is the distance of the donor charge from the interface at z=-d (donor charge is assimed to be at (0,0,0)
epsilon is the dielectric constant of the host material where the donor charge resides
epsilonrel is (e2-e1)/(e2+e1) which quantifies the correction. A barrier of height vbarr is placed at z=-d to confine the particle near but inside the interface.
To avoid sugularities the inteface os shifted by -del away from -d. 

 
Ian Galbraith , Heriot Watt University, Scotland. I.Galbraith@hw.ac.uk
*/

PrimImageAction::PrimImageAction(const double d, const double del, const double epsilon , const double epsilonrel, const double vbarr,
  const SimulationInfo &simInfo, int ndim, const Species &species)
: tau(simInfo.getTau()), d(d), del(del),  epsilon(epsilon), epsilonrel(epsilonrel), vbarr(vbarr), ndim(ndim),
  ifirst(species.ifirst), npart(species.count){
}

double PrimImageAction::getActionDifference(const SectionSamplerInterface& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  double ktstride = tau*nStride;
  double rho2=0; 
	double vlocal=0;
  double z=0;
	double denom1=0;
	double denom2=0;


  for (int islice=nStride; islice<nSlice-nStride; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
		if (i<ifirst || i>=ifirst+npart) continue;
      // Add action for moving beads.
      Vec delta=movingBeads(iMoving,islice);
      cell.pbc(delta);
	
	    rho2=delta[0]*delta[0] +  delta[1]*delta[1];
		z=delta[2];
		if (z < -d) {
			vlocal= vbarr;
		}else {
			vlocal=0.0;
		}

		denom1=sqrt(rho2+(z+2*(d+del))*(z+2*(d+del)));
		denom2= (z+d+del);
		
		
		deltaAction+=( vlocal + ( 2/denom1 - 0.5/denom2 )*epsilonrel/epsilon )*ktstride;
      // Subtract action for old beads.
      delta=sectionBeads(i,islice);
      cell.pbc(delta);
		rho2=delta[0]*delta[0] +  delta[1]*delta[1];
		z=delta[2];
		if (z < -d) {
			vlocal= vbarr;
		}else {
			vlocal=0.0;
		}

		denom1=sqrt(rho2+(z+2*(d+del))*(z+2*(d+del)));
		denom2= (z+d+del);
		
		
		deltaAction-= ( vlocal + ( 2/denom1 - 0.5/denom2 )*epsilonrel/epsilon ) *ktstride;	
    }
  }
  return deltaAction;
}

double PrimImageAction::getTotalAction(const Paths& paths, 
    const int level) const {
  return 0;
}

void PrimImageAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double& u, double& utau, double& ulambda, Vec &fm, Vec &fp) const {
	return;
}
