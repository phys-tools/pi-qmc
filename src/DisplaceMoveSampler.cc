// $Id: MultiLevelSampler.cc 22 2009-5-18 Saad khairallah$
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
#include "DisplaceMoveSampler.h"
#include "stats/AccRejEstimator.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "Action.h"
#include "DoubleAction.h"
#include "Mover.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include "ParallelPaths.h"
#include "Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include <sstream>
#include <string>


DisplaceMoveSampler::DisplaceMoveSampler(const int nmoving, Paths& paths, const double dist, const double freq,
					 ParticleChooser& particleChooser, Mover& mover, Action* action, 
					 const int nrepeat, const BeadFactory& beadFactory,
					 const MPIManager* mpi)
  : nmoving(nmoving), dist(dist), mover(mover), cell(paths.getSuperCell()), action(action), freq(freq),
    movingIndex(new IArray(nmoving)), identityIndex(nmoving),
    particleChooser(particleChooser), paths(paths), accRejEst(0), nrepeat(nrepeat), mpi(mpi), 
    beadFactory(beadFactory){
 
  nslice = paths.getnprocSlice();
  pathsBeads = beadFactory.getNewBeads(paths.getNPart(), nslice);
  for (int i=0; i<nmoving; ++i) identityIndex(i)=i;
}

DisplaceMoveSampler::~DisplaceMoveSampler() {
  delete movingBeads;
  delete movingIndex;
  delete &particleChooser;
}



void DisplaceMoveSampler::run() {
  
  // run with a probability equal to freq 
#ifdef ENABLE_MPI
  if ( mpi && (mpi->getNWorker()) > 1) {
    double rd = RandomNumGenerator::getRand() ;
    mpi->getWorkerComm().Bcast(&rd, 1, MPI::DOUBLE, 0); 
    if (rd > freq) return;
  } else {
    if (RandomNumGenerator::getRand()>freq) return ;
  }
#else
  if (RandomNumGenerator::getRand()>freq) return ;
#endif
  
  const Permutation &pathsPermutation(paths.getPermutation());
  
  
  // Select particles that are not permuting to move
  iFirstSlice = paths.getLowestSampleSlice(0,false);
  paths.getBeads(iFirstSlice,*pathsBeads);
  
  action->initialize(*this);
  
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    movingIndex->resize(nmoving);  
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
      identityIndex.resize(imovingNonPerm);
      
      // Copy old coordinate to the moving coordinate
      movingBeads = beadFactory.getNewBeads(imovingNonPerm, nslice);
      for (int i=0; i<imovingNonPerm; ++i) identityIndex(i)=i;
      for (int islice=0; islice<nslice; ++islice) { 
	pathsBeads->copySlice(*movingIndex,islice,*movingBeads,identityIndex,islice);
      }
      
      if (tryMove(imovingNonPerm)) continue;
      delete movingBeads;
    }
  }
  
}




bool DisplaceMoveSampler::tryMove(int imovingNonPerm) {
 
  accRejEst->tryingMove(0);
  double l = mover.makeMove(*this);

  // Evaluate the change in action.
  double deltaAction=(action==0)?0:action->getActionDifference(*this, imovingNonPerm);

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

  // Put moved beads in paths beads.
  for (int islice=0; islice<nslice; ++islice) {
    movingBeads->copySlice(identityIndex,islice,
			   *pathsBeads,*movingIndex,islice);
  }

  
  Permutation perm(paths.getNPart()); 
  paths.putBeads(iFirstSlice,*pathsBeads,perm);

  return true;
}



void DisplaceMoveSampler::setAction(Action* act, const int level) {action=act;}



AccRejEstimator* 
DisplaceMoveSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": moving " << nmoving
           << " " << particleChooser.getName();
  return accRejEst=new AccRejEstimator(longName.str(),1);
}
