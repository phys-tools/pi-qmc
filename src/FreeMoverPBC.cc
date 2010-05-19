// $Id: FreeMoverPBC.cc 250 2010-04-06 23:12:06Z saadAK $
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
#include "FreeMoverPBC.h"
#include "Beads.h"
#include "MultiLevelSampler.h"
#include "RandomNumGenerator.h"
#include <blitz/tinyvec.h>
#include "SuperCell.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"
#include <cmath>

FreeMoverPBC::FreeMoverPBC(const double lam, const int npart, const double tau)
  : lambda(npart), tau(tau) {
  lambda=lam; 
}

FreeMoverPBC::FreeMoverPBC(const SimulationInfo& simInfo, const int maxlevel,
  const double pgDelta)
  : lambda(simInfo.getNPart()), tau(simInfo.getTau()),
    pg(maxlevel+1,simInfo.getNSpecies(),NDIM),
    length(simInfo.getSuperCell()->a), specIndex(simInfo.getNPart()) {
  for (int i=0; i<simInfo.getNPart(); ++i) {
    const Species* species=&simInfo.getPartSpecies(i);
    lambda(i)=0.5/species->mass;
    for (int j=0; j<simInfo.getNSpecies(); ++j) {
      if (species==&simInfo.getSpecies(j)) {specIndex(i)=j; break;}
    }
  }
  const int nspec=simInfo.getNSpecies();
  pg.resize(maxlevel+1,nspec,NDIM);
  maxNImage = 27;
  for (int ilevel=0; ilevel<=maxlevel; ++ilevel) {
    for (int ispec=0; ispec<nspec; ++ispec) {
      double alpha=simInfo.getSpecies(ispec).mass/(tau*pow(2,ilevel));
      int nimage = 1;
      for (int idim=0; idim<NDIM; ++idim) {
        double l = length(idim);
        if (l*l*alpha<60) {
          pg(ilevel,ispec,idim) 
              = new PeriodicGaussian(alpha,l,(int)(100*l*sqrt(alpha)));
          nimage *= 1+2*ceil(5./l*sqrt(alpha));
        } else {
          pg(ilevel,ispec,idim)=0;
        }
      }
      if (nimage>maxNImage) maxNImage=nimage;
    }
  }
  //std::cout << "maxNImage=" << maxNImage << std::endl;
  weight.resize(maxNImage);
}

FreeMoverPBC::~FreeMoverPBC() {
  for (PGArray::iterator i=pg.begin();  i!=pg.end(); ++i) delete *i;
}

double FreeMoverPBC::makeMove(MultiLevelSampler& sampler, const int level) {
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
      // First sample periodic image to determine midpoint .
      IVec nimage;
      for (int idim=0; idim<NDIM; ++idim) {
        nimage(idim) = 1+2*ceil(5.*sigma/length(idim));
      } 
      //std::cout << "level, nimage" << level << ", " << nimage << std::endl;
      Vec midpoint=movingBeads.delta(iMoving,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      midpoint+=movingBeads(iMoving,islice-nStride);
      // Now sample about midpoint
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
double FreeMoverPBC::makeDelayedMove(MultiLevelSampler& sampler, const int level) {
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
      Vec delta=rejectedBeads(iMoving,islice); delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= (*pg(level,ispec,idim))(fabs(delta[idim]));
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
      delta = movingBeads(iMoving,islice);delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        if (pg(level,ispec,idim)) {
          temp *= (*pg(level,ispec,idim))(fabs(delta[idim]));
        } else {
	  toldOverTnew -=  delta[idim]*delta[idim]*inv2Sigma2;
        }
      }
      toldOverTnew+=log(temp);
    }
  }
  return toldOverTnew; //Return the log of the probability.
}
