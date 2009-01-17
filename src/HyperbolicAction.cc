// $Id: HyperbolicAction.cc,v 1.8 2007/10/03 12:53:56 jshumwa Exp $
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
#include "HyperbolicAction.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"
#include <gsl/gsl_sf_bessel.h>

HyperbolicAction::HyperbolicAction(const SimulationInfo& simInfo, 
                                   const int maxlevel) 
  : tau(simInfo.getTau()) {
  //const int nspec=simInfo.getNSpecies();
  /*for (int i=0; i<simInfo.getNPart(); ++i) {
    const Species* species=&simInfo.getPartSpecies(i);
    lambda(i)=0.5/species->mass;
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
        if (l*l*alpha<10) {
          pg(ilevel,ispec,idim)=0;//new PeriodicGaussian(alpha,l,
                                   //            (int)((0.5*l)/deltaPG+0.01));
        } else {
          pg(ilevel,ispec,idim)=0;
        }
      }
    }
  }
 */
}

HyperbolicAction::~HyperbolicAction() {
}

double HyperbolicAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  double teff=tau*nStride; 
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  const double mass=0.067;
  const double alpha=27.211396/1.525;
  double deltaAction=0;
  for (int islice=nStride; islice<nSlice; islice+=nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      // Add action for moving beads.
      Vec delta=movingBeads.delta(iMoving,islice,-nStride);
      cell.pbc(delta);
      double temp = 1.0;
      double z=teff/(2*alpha)*sqrt(1+(2*mass*alpha*dot(delta,delta))
                                    /(teff*teff));
      temp*=gsl_sf_bessel_Kn(2,z)/(z*z);
      temp = 1.0/temp;
      // Subtract action for old beads.
      delta=sectionBeads.delta(i,islice,-nStride);
      cell.pbc(delta);
      z=teff/(2*alpha)*sqrt(1+(2*mass*alpha*dot(delta,delta))/(teff*teff));
      temp*=gsl_sf_bessel_Kn(2,z)/(z*z);
      deltaAction+=log(temp);
    }
  }
  return deltaAction;
}

double HyperbolicAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void HyperbolicAction::getBeadAction(const Paths& paths, int ipart, int islice,
     double &u, double &utau, double &ulambda, Vec &fm, Vec &fp) const {
  u=utau=ulambda=0; fm=0.; fp=0.;
  const double mass=0.067;
  const double alpha=27.211396/1.525;
  Vec delta = paths.delta(ipart,islice,-1);
  double z=tau/(2*alpha)*sqrt(1+(2*mass*alpha*dot(delta,delta))/(tau*tau));
  utau = -1./(2.*alpha)-1./tau+tau/(2*alpha*alpha*z*z)
         +(gsl_sf_bessel_Kn(1,z)+gsl_sf_bessel_Kn(3,z))*tau
         /(gsl_sf_bessel_Kn(2,z)*8*alpha*alpha*z);
}
