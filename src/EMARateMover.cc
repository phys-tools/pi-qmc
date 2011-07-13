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

EMARateMover::EMARateMover(const SimulationInfo& simInfo, const int maxlevel,
  const double pgDelta)
  : PermutationChooser(2), ParticleChooser(2),
    lambda(simInfo.getNPart()), tau(simInfo.getTau()) {
  for (int i=0; i<simInfo.getNPart(); ++i) {
    const Species* species=&simInfo.getPartSpecies(i);
    lambda(i)=0.5/species->mass;
  }
}

EMARateMover::~EMARateMover() {
}

double EMARateMover::makeMove(MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const blitz::Array<int,1>& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  blitz::Array<Vec,1> gaussRand(nMoving); gaussRand=0.0;
  double toldOverTnew=0;  
  forwardProb = 0;
  // Make radiating vs. direct choice at highest level.
  std::cout << nStride << ", " << nSlice << std::endl;
  if (nStride+1==nSlice) {
    std::cout << "Choosing radiating vs. diagonal." << std::endl; 
    double i = index(0);
    double j = index(1);
    Vec re1 = movingBeads(0,0);
    Vec re2 = movingBeads(0,nSlice-1);
    Vec rh1 = movingBeads(1,0);
    Vec rh2 = movingBeads(1,nSlice-1);
    Vec deltae = re2-re1; cell.pbc(deltae);
    Vec deltah = rh2-rh1; cell.pbc(deltah);
    Vec delta1 = re1-rh1; cell.pbc(delta1);
    Vec delta2 = re2-rh2; cell.pbc(delta2);
    double inv2sigma2e = 0.5/(lambda(i)*tau*nStride); 
    double inv2sigma2h = 0.5/(lambda(j)*tau*nStride); 
    double inv2sigma21 = 1./((lambda(i)+lambda(j))*tau*nStride); 
    double te = dot(deltae,deltae)*inv2sigma2e;
    double th = dot(deltah,deltah)*inv2sigma2h;
    double t1 = dot(delta1,delta1)*inv2sigma21;
    double t2 = dot(delta2,delta2)*inv2sigma21;
    std::cout << t1 << ", " << t2 << ", " << te << ", " << th << std::endl;
    double pRad = exp(-(t1+t2));
    double pDiag = exp(-(te+te));
    std::cout << pRad << ", " << pDiag << std::endl;
    if (RandomNumGenerator::getRand()>pDiag/(pDiag+pRad)) {
      std::cout << "Trying radiating move." << std::endl;
    } else {
      std::cout << "Trying diagonal move." << std::endl;
    }
  } {

  for (int islice=nStride; islice<nSlice-nStride; islice+=2*nStride) {
    RandomNumGenerator::makeGaussRand(gaussRand);
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      double sigma = sqrt(lambda(i)*tau*nStride); 
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
      for (int idim=0;idim<NDIM;++idim) {
        double tmp = delta[idim]*delta[idim]*inv2Sigma2;
        forwardProb += tmp;
        toldOverTnew += tmp;
      }

      // Calculate and add reverse transition probability.
      midpoint=sectionBeads.delta(i,islice+nStride,-2*nStride);
      cell.pbc(midpoint)*=0.5;
      midpoint+=sectionBeads(i,islice-nStride);
      delta=sectionBeads(i,islice); delta-=midpoint;
      cell.pbc(delta);
      for (int idim=0;idim<NDIM;++idim) {
        toldOverTnew -= delta[idim]*delta[idim]*inv2Sigma2;
      }
    }
  }
  }
  return toldOverTnew; //Return the log of the probability.
}


// Delayed Rejection 
double EMARateMover::makeDelayedMove(MultiLevelSampler& sampler,
    const int level) {
  return 1.0;
}
