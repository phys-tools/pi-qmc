// $Id$
/*  Copyright (C) 2010 John B. Shumway, Jr. and Saad A. Khairallah

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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "EMARateMover.h"
#include "Beads.h"
#include "MultiLevelSampler.h"
#include "RandomNumGenerator.h"
#include <blitz/tinyvec.h>
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"
#include <cmath>

EMARateMover::EMARateMover(const double lam, const int npart, const double tau)
  : lambda(npart), tau(tau) {
  lambda=lam; 
}

EMARateMover::EMARateMover(const SimulationInfo& simInfo, const int maxlevel,
  const double pgDelta)
  : lambda(simInfo.getNPart()), tau(simInfo.getTau()),
    pg(maxlevel+1,simInfo.getNSpecies(),NDIM), specIndex(simInfo.getNPart()) {
  for (int i=0; i<simInfo.getNPart(); ++i) {
    const Species* species=&simInfo.getPartSpecies(i);
    lambda(i)=0.5/species->mass;
    for (int j=0; j<simInfo.getNSpecies(); ++j) {
      if (species==&simInfo.getSpecies(j)) {specIndex(i)=j; break;}
    }
  }
  const int nspec=simInfo.getNSpecies();
  pg.resize(maxlevel+1,nspec,NDIM);
  for (int ilevel=0; ilevel<=maxlevel; ++ilevel) {
    for (int ispec=0; ispec<nspec; ++ispec) {
      double alpha=simInfo.getSpecies(ispec).mass/(tau*pow(2,ilevel));
      for (int idim=0; idim<NDIM; ++idim) {
        double l=(*simInfo.getSuperCell())[idim];
        if (l*l*alpha<60) {
          pg(ilevel,ispec,idim)=new PeriodicGaussian(alpha,l,
                                                     (int)(100*l*sqrt(alpha)));
        } else {
          pg(ilevel,ispec,idim)=0;
        }
      }
    }
  }
}

EMARateMover::~EMARateMover() {
  for (PGArray::iterator i=pg.begin();  i!=pg.end(); ++i) delete *i;
}

double EMARateMover::makeMove(MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  double factor = sampler.getFactor(); 
  const blitz::Array<int,1>& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  blitz::Array<Vec,1> gaussRand(nMoving); gaussRand=0.0;
  double toldOverTnew=0;  
  forwardProb = 0;
  for (int islice=nStride; islice<nSlice-nStride; islice+=2*nStride) {
    RandomNumGenerator::makeGaussRand(gaussRand);
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      const int ispec=specIndex(i);
      double sigma = factor*sqrt(lambda(i)*tau*nStride); 
      double inv2Sigma2 = 0.5/(sigma*sigma);
      // Calculate the new position.
      Vec midpoint=movingBeads.delta(iMoving,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      midpoint+=movingBeads(iMoving,islice-nStride);
      Vec delta = gaussRand(iMoving); delta*=sigma;
      (movingBeads(iMoving,islice)=midpoint)+=delta;
      cell.pbc(movingBeads(iMoving,islice));
      // Calculate 1/transition probability for move 
      cell.pbc(delta);
      double temp = 1.0;
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= (*pg(level,ispec,idim))(fabs(delta[idim]));
        } else {
	  double tmp = delta[idim]*delta[idim]*inv2Sigma2;
	  forwardProb += tmp;
	  toldOverTnew += tmp;
	}
      }
      temp = 1.0/temp;    
      forwardProb +=log(temp);

      // Calculate and add reverse transition probability.
      midpoint=sectionBeads.delta(i,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      midpoint+=sectionBeads(i,islice-nStride);
      delta=sectionBeads(i,islice); delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= (*pg(level,ispec,idim))(fabs(delta[idim]));
        } else {
          toldOverTnew -= delta[idim]*delta[idim]*inv2Sigma2;
        }
      }
      toldOverTnew+=log(temp);
    }
  }
  return toldOverTnew; //Return the log of the probability.
}


// Delayed Rejection 
double EMARateMover::makeDelayedMove(MultiLevelSampler& sampler,
    const int level) {
  return 1.0;
}
