// $Id: DoubleDisplaceSampler.cc 22 2009-5-18 Saad khairallah$
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
#include "stats/MPIManager.h"
#include "DoubleDisplaceMoveSampler.h"
#include "stats/AccRejEstimator.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "Action.h"
#include "DoubleAction.h"
#include "Mover.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include "Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include <sstream>
#include <string>


DoubleDisplaceMoveSampler::DoubleDisplaceMoveSampler(const int nmoving, Paths& paths, const double dist, 
						     const double freq, ParticleChooser& particleChooser, 
						     Mover& mover, Action* action, const int nrepeat, 
						     const BeadFactory& beadFactory, const MPIManager* mpi, 
						     DoubleAction* doubleAction)
  : DisplaceMoveSampler(nmoving, paths, dist, freq, particleChooser, mover, action, nrepeat, beadFactory, mpi),doubleAction(doubleAction){
  movingIndex2 = new IArray(nmoving);
  nslice = paths.getnprocSlice();  if (mpi->getNWorker() ==1 && mpi->getNClone()!=0){ nslice/=2;  pathsBeads = beadFactory.getNewBeads(paths.getNPart(), nslice);}
  pathsBeads2 = beadFactory.getNewBeads(paths.getNPart(), nslice); 
  pathsBeads1 = beadFactory.getNewBeads(paths.getNPart(), nslice);
}

DoubleDisplaceMoveSampler::~DoubleDisplaceMoveSampler() {
  delete movingBeads;  delete movingBeads2; delete movingBeads1;
  delete movingIndex;  delete movingIndex2;
  delete &particleChooser;
}



void DoubleDisplaceMoveSampler::run() {
 
  // Run with a probability equal to freq 
#ifdef ENABLE_MPI
  if ( mpi && (mpi->getNWorker()) > 1) {
    double rd = RandomNumGenerator::getRand() ;
    mpi->getWorkerComm().Bcast(&rd, 1, MPI::DOUBLE, 0); 
    if (rd > freq || nrepeat ==0) return;
  } else {
    if (RandomNumGenerator::getRand()>freq || nrepeat ==0) return ;
  }
#else
  if (RandomNumGenerator::getRand()>freq || nrepeat ==0) return ;
#endif

 
  const Permutation &pathsPermutation(paths.getPermutation());

  iFirstSlice = paths.getLowestSampleSlice(0,false);
  paths.getBeads(iFirstSlice,*pathsBeads1);

  iFirstSlice2 = (iFirstSlice + paths.getNSlice()/2)%paths.getNSlice();
  paths.getBeads(iFirstSlice2,*pathsBeads2);
  
  action->initialize(*this);
  doubleAction->initialize(*this);

  // Select particles that are not permuting to move and make nrepeat displace moves
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    movingIndex->resize(nmoving); 
    movingIndex2->resize(nmoving); 
    identityIndex.resize(nmoving);  
    
    particleChooser.chooseParticles();
    int imovingNonPerm = 0;         
    for (int i=0; i<nmoving; ++i) {
      int j = particleChooser[i];
      int jperm = pathsPermutation[j];
      
      if (j==jperm){
	(*movingIndex)(imovingNonPerm)=j;
	imovingNonPerm++;
      }
    } 
    
 
#ifdef ENABLE_MPI
    if ( mpi && (mpi->getNWorker()) > 1) {
      mpi->getWorkerComm().Bcast(&(*movingIndex)(0), nmoving, MPI::INT, 0); 
      mpi->getWorkerComm().Bcast(&imovingNonPerm, 1, MPI::INT, 0); 
    }
#endif
    
    if (imovingNonPerm !=0 ) {
      movingIndex->resizeAndPreserve(imovingNonPerm);      
      movingIndex2->resizeAndPreserve(imovingNonPerm);
      identityIndex.resize(imovingNonPerm);
      
      // Copy old coordinate to the moving coordinate
      for (int i=0; i<imovingNonPerm; ++i) identityIndex(i)=i;
      
      movingBeads1 = beadFactory.getNewBeads(imovingNonPerm, nslice);
      movingBeads2 = beadFactory.getNewBeads(imovingNonPerm, nslice);

      ////
      for (int i=0; i<imovingNonPerm; i++) (*movingIndex2)(i) = (*movingIndex)(i);//////////<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ////
      for (int islice=0; islice<nslice; ++islice) { 
	pathsBeads1->copySlice(*movingIndex,islice,*movingBeads1,identityIndex,islice);
	pathsBeads2->copySlice(*movingIndex2,islice,*movingBeads2,identityIndex,islice); /////<<<<<<<<<<<<<<<<<<<<<<<
      }
      
      if (tryMove(imovingNonPerm)) continue;
   
      delete movingBeads2; 
      delete movingBeads1;
    }
  }

}




bool DoubleDisplaceMoveSampler::tryMove(int imovingNonPerm) {

  //suggest a move
  accRejEst->tryingMove(0);
  double l = mover.makeMove(*this);

  // Evaluate the change in action. 
  activateSection(2);
  double deltaAction=(action==0)?0:action->getActionDifference(*this, imovingNonPerm);
  activateSection(1);
  deltaAction+=(action==0)?0:action->getActionDifference(*this, imovingNonPerm);
  deltaAction+=(action==0)?0:doubleAction->getActionDifference(*this, imovingNonPerm);

#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    double netDeltaAction=0;
    mpi->getWorkerComm().Reduce(&deltaAction,&netDeltaAction,1,MPI::DOUBLE,MPI::SUM,0);
    double acceptProb = exp(-netDeltaAction);
    bool acceptReject = RandomNumGenerator::getRand()>acceptProb;
    mpi->getWorkerComm().Bcast(&acceptReject,1,MPI::CHAR,0); 
    if (acceptReject) return false;
  }else
    {
      double acceptProb=exp(-deltaAction); 
      if (RandomNumGenerator::getRand()>acceptProb) return false;
    }
#else
  double acceptProb=exp(-deltaAction);
  if (RandomNumGenerator::getRand()>acceptProb) return false;
#endif
  
  accRejEst->moveAccepted(0);
  
  // Move accepted.
  action->acceptLastMove();  
  doubleAction->acceptLastMove();

  // Put moved beads in paths beads.
  for (int islice=0; islice<nslice; ++islice) {
    movingBeads1->copySlice(identityIndex,islice,*pathsBeads1,*movingIndex,islice);
    movingBeads2->copySlice(identityIndex,islice,
			    *pathsBeads2,*movingIndex2,islice);///////////<<<<<<<<<<<<<<<<<<<<<<<<<<,
  } 

  // Append the currenrmutation to section permutation.
  Permutation perm1(paths.getNPart()); 
  Permutation perm2(paths.getNPart()); 

  paths.putDoubleBeads(iFirstSlice,*pathsBeads1,perm1,
		       iFirstSlice2,*pathsBeads2,perm2);
  
  return true;
}

void DoubleDisplaceMoveSampler::activateSection(const int i) {
  if (i==1) {
    pathsBeads = pathsBeads1;
    movingBeads = movingBeads1;
  } else {
    pathsBeads = pathsBeads2;
    movingBeads = movingBeads2;
  }
}

void DoubleDisplaceMoveSampler::setAction(Action* act, const int level) {action=act;}



AccRejEstimator* 
DoubleDisplaceMoveSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
 
  longName << name << ": moving " << nmoving
           << " " << particleChooser.getName();
  return accRejEst=new AccRejEstimator(longName.str(),1);
}
