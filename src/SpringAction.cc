// $Id$
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
#include "SpringAction.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "DisplaceMoveSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"

SpringAction::SpringAction(const SimulationInfo& simInfo, const int maxlevel,
  const double deltaPG) 
  : lambda(simInfo.getNPart()), tau(simInfo.getTau()), 
    pg(maxlevel+1,simInfo.getNSpecies(),NDIM),
    specIndex(simInfo.getNPart()), isStatic(simInfo.getNPart()) {
  const int nspec=simInfo.getNSpecies();
  for (int i=0; i<simInfo.getNPart(); ++i) {
    const Species* species=&simInfo.getPartSpecies(i);
    lambda(i)=0.5/species->mass;
    isStatic(i)=species->isStatic;
    for (int j=0; j<nspec; ++j) {
      if (species==&simInfo.getSpecies(j)) {specIndex(i)=j ;break;}
    }
  }
  pg.resize(maxlevel+1,nspec,NDIM);
  for (int ilevel=0; ilevel<=maxlevel; ++ilevel) {
    for (int ispec=0; ispec<nspec; ++ispec) {
      double alpha=simInfo.getSpecies(ispec).mass/(2*tau*pow(2,ilevel));
      for (int idim=0; idim<NDIM; ++idim) {
        double l=(*simInfo.getSuperCell())[idim];
        if (l*l*alpha<60) {
          int ngridpts = (int)(100*l*sqrt(alpha));
          pg(ilevel,ispec,idim) = new PeriodicGaussian(alpha,l,ngridpts);
          if (ilevel==0) {
            std::cout << "WARNING: Periodic Gaussian will give incorrect"
              << " energy.\n         Decrease tau or fix SpringAction."
              << std::endl;
          }
        } else {
          pg(ilevel,ispec,idim)=0;
        }
      }
    }
  }
}

SpringAction::~SpringAction() {
  for (PGArray::iterator i=pg.begin();  i!=pg.end(); ++i) delete *i;
}

double SpringAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      if (isStatic(i)) continue;
      const int ispec=specIndex(i);
      const double inv2Sigma2 = 0.25/(lambda(i)*tau*nStride);
      // Add action for moving beads.
      Vec delta=movingBeads.delta(iMoving,islice,-nStride);
      cell.pbc(delta);
      double temp = 1.0;
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= (*pg(level,ispec,idim))(fabs(delta[idim]));
        } else {
          deltaAction+=delta[idim]*delta[idim]*inv2Sigma2;
        }
      }
      temp = 1.0/temp;
      // Subtract action for old beads.
      delta=sectionBeads.delta(i,islice,-nStride);
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= (*pg(level,ispec,idim))(fabs(delta[idim]));
        } else {
          deltaAction-=delta[idim]*delta[idim]*inv2Sigma2;
        }
      }
      deltaAction+=log(temp);
    }
  } 
  return deltaAction;
}


double SpringAction::getActionDifference(const Paths &paths, 
   const VArray &displacement, int nmoving, const IArray &movingIndex, 
   int iFirstSlice, int iLastSlice) {
 return 0; //No change in action for uniform displacements of particles.
}


double SpringAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void SpringAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double &u, double &utau, double &ulambda, Vec &fm, Vec &fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  if (isStatic(ipart)) return;
  Vec delta = paths.delta(ipart,islice,-1);
  fm-=(delta/(2*lambda(ipart)*tau));
  utau = NDIM/(2*tau) - dot(delta,delta)/(4.0*lambda(ipart)*tau*tau);   

  delta = paths.delta(ipart,islice,1);
  fp-=(delta/(2*lambda(ipart)*tau));

}
