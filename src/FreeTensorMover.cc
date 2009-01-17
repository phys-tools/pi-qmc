// $Id: FreeTensorMover.cc,v 1.7 2008/11/25 13:11:21 jshumwa Exp $
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
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "FreeTensorMover.h"
#include "Beads.h"
#include "MultiLevelSampler.h"
#include "RandomNumGenerator.h"
#include <blitz/tinyvec.h>
#include "SuperCell.h"
#include "SimulationInfo.h"
#include <cmath>

FreeTensorMover::FreeTensorMover(const SimulationInfo& simInfo)
  : lambda(simInfo.getNPart()), tau(simInfo.getTau()) {
  for (int i=0; i<simInfo.getNPart(); ++i) {
    lambda(i)=0.5;
    lambda(i)/=(*simInfo.getPartSpecies(i).anMass);
  }
}

double FreeTensorMover::makeMove(MultiLevelSampler& sampler, const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const blitz::Array<int,1>& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  blitz::Array<Vec,1> gaussRand(nMoving);
  double toldOverTnew=0;
  for (int islice=nStride; islice<nSlice-nStride; islice+=2*nStride) {
    RandomNumGenerator::makeGaussRand(gaussRand);
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      const int i=index(iMoving);
      Vec& lambdai=lambda(i);
      double sigmax=sqrt(lambdai[0]*tau*nStride);
      double sigmay=sqrt(lambdai[1]*tau*nStride);
      double sigmaz=sqrt(lambdai[2]*tau*nStride);
      double inv2Sigma2x=0.5/(sigmax*sigmax);
      double inv2Sigma2y=0.5/(sigmay*sigmay);
      double inv2Sigma2z=0.5/(sigmaz*sigmaz);
      // Calculate the new position.
      Vec midpoint=movingBeads(iMoving,islice+nStride);
      cell.pbc(midpoint-=movingBeads(iMoving,islice-nStride))*=0.5;
      midpoint+=movingBeads(iMoving,islice-nStride);
      Vec delta = gaussRand(iMoving);
      delta[0]*=sigmax; delta[1]*=sigmay; delta[2]*=sigmaz;
      (movingBeads(iMoving,islice)=midpoint)+=delta;
      cell.pbc(movingBeads(iMoving,islice));
      // Add transition probability for move.
      toldOverTnew+=delta[0]*delta[0]*inv2Sigma2x
                   +delta[1]*delta[1]*inv2Sigma2y
                   +delta[2]*delta[2]*inv2Sigma2z;
      // Calculate and add reverse transition probability.
      midpoint=sectionBeads(i,islice+nStride);
      cell.pbc(midpoint-=sectionBeads(i,islice-nStride))*=0.5;
      midpoint+=sectionBeads(i,islice-nStride);
      delta=sectionBeads(i,islice); delta-=midpoint;
      cell.pbc(delta);
      toldOverTnew-=delta[0]*delta[0]*inv2Sigma2x
                   +delta[1]*delta[1]*inv2Sigma2y
                   +delta[2]*delta[2]*inv2Sigma2z;
    }
  }
  ///toldOverTnew=exp(toldOverTnew);
  return toldOverTnew;  // return log of probability.
}
