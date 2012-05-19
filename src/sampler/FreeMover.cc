// $Id$
/*  Copyright (C) 2004-2006 John B. Shumway, Jr. and Saad A. Khairallah

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
#include "FreeMover.h"
#include "Beads.h"
#include "MultiLevelSampler.h"
#include "util/RandomNumGenerator.h"
#include <cstdlib>
#include <blitz/tinyvec.h>
#include "util/SuperCell.h"
#include "SimulationInfo.h"
#include "util/PeriodicGaussian.h"
#include <cmath>

FreeMover::FreeMover(const double lam, const int npart, const double tau)
  : lambda(npart), tau(tau) {
  lambda=lam; 
}

FreeMover::FreeMover(const SimulationInfo& simInfo, const int maxlevel,
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
        double length = (*simInfo.getSuperCell())[idim];
        if (PeriodicGaussian::numberOfTerms(alpha, length) < 16) {
          pg(ilevel,ispec,idim) = new PeriodicGaussian(alpha,length);
        } else {
          pg(ilevel,ispec,idim) = 0;
        }
      }
    }
  }
}

FreeMover::~FreeMover() {
  for (PGArray::iterator i=pg.begin();  i!=pg.end(); ++i) delete *i;
}

double FreeMover::makeMove(MultiLevelSampler& sampler, const int level) {
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
      cell.pbc(midpoint);
      Vec delta = gaussRand(iMoving); 
      delta*=sigma;  
      cell.pbc(delta);
      (movingBeads(iMoving,islice)=midpoint)+=delta;
      cell.pbc(movingBeads(iMoving,islice));
      // Calculate 1/transition probability for move 
      double temp = 1.0;
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
            temp *= pg(level, ispec, idim)->evaluate(delta[idim]);
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
      cell.pbc(midpoint);
      delta=sectionBeads(i,islice); 
      delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= pg(level,ispec,idim)->evaluate(delta[idim]);
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
double FreeMover::makeDelayedMove(MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& rejectedBeads=sampler.getRejectedBeads();
  Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=movingBeads.getNSlice(); 
  double factor = sampler.getFactor(); 
  const blitz::Array<int,1>& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();

  double toldOverTnew=0; 
  forwardProb = 0;
  for (int islice=nStride; islice<nSlice-nStride; islice+=2*nStride) {
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      const int ispec=specIndex(i);
      double sigma = factor*sqrt(lambda(i)*tau*nStride); 
      double inv2Sigma2 = 0.5/(sigma*sigma);
      double temp = 1.0;
      // Calculate forward transition probability C2B 
      Vec midpoint=movingBeads.delta(iMoving,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      midpoint+=movingBeads(iMoving,islice-nStride);
      cell.pbc(midpoint);
      Vec delta=rejectedBeads(iMoving,islice); 
      delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= pg(level,ispec,idim)->evaluate(delta[idim]);
        } else {
	  double tmp = delta[idim]*delta[idim]*inv2Sigma2;
	  forwardProb -= tmp;
          toldOverTnew += tmp;
        }
      }

      forwardProb +=log(temp);         
      temp = 1.0/temp;

      // Calculate inverse transition probability B2C i.e. T(R1->R2)=numerator
      midpoint=rejectedBeads.delta(iMoving,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      midpoint+=rejectedBeads(iMoving,islice-nStride);
      cell.pbc(midpoint);
      delta = movingBeads(iMoving,islice);
      delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= pg(level,ispec,idim)->evaluate(delta[idim]);
        } else {
	  toldOverTnew -=  delta[idim]*delta[idim]*inv2Sigma2;
        }
      }
      toldOverTnew+=log(temp);
    }
  }
  return toldOverTnew; //Return the log of the probability.
}
