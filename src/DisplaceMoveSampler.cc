// $Id: MultiLevelSampler.cc 22 2009-5-18 Saad khairallah$
/*  Copyright (C) 2009 John B. Shumway, Jr.

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
#include "Action.h"
#include "UniformMover.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include "Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include "SuperCell.h"
#include <sstream>
#include <string>

DisplaceMoveSampler::DisplaceMoveSampler(int nmoving, int nrepeat,
  Paths& paths, ParticleChooser& particleChooser, const UniformMover& mover, 
  Action* action, const MPIManager* mpi)
  : nmoving(nmoving), nrepeat(nrepeat),
    nslice(paths.getNUniqueSlice()),
    iFirstSlice(paths.getLowestSampleSlice(1,true)),
    displacement(nmoving), movingIndex(nmoving), 
    particleChooser(particleChooser), mover(mover), 
    paths(paths), cell(paths.getSuperCell()), action(action),
    accRejEst(0),  mpi(mpi) {
}

DisplaceMoveSampler::~DisplaceMoveSampler() {
  delete &particleChooser;
}


void DisplaceMoveSampler::run() {
  int workerID = (mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  int nworker = (mpi)?mpi->getNWorker():1;
#endif
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    
    // Select particles at random and reject move if any of the selected
    // particles are permuting.
    bool reject=false;
    if (workerID==0) {
      particleChooser.chooseParticles();
      Permutation permutation=paths.getPermutation();
      for (int i=0; i<nmoving; ++i) {
        int j = particleChooser[i];
        if ( j == permutation[j] ){
          movingIndex(i)=j;
        } else {
          reject=true;
        }
      } 
    }
#ifdef ENABLE_MPI
    if (nworker > 1) {
      mpi->getWorkerComm().Bcast(&reject, 1, MPI::CHAR, 0); 
    }
#endif 
    if (reject) continue;
#ifdef ENABLE_MPI
    if ( nworker > 1) {
      mpi->getWorkerComm().Bcast(&(movingIndex)(0), nmoving, MPI::INT, 0); 
    }
#endif
    tryMove();
  }
}


bool DisplaceMoveSampler::tryMove() {
 
  accRejEst->tryingMove(0);
  mover.makeMove(displacement,nmoving);

  // Evaluate the change in action.
  double deltaAction = (action==0) ?  0 
   : action->getActionDifference(paths,displacement,nmoving,movingIndex,
                                 iFirstSlice,nslice);

#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    double netDeltaAction=0;
    mpi->getWorkerComm().Reduce(&deltaAction,&netDeltaAction,
                                1,MPI::DOUBLE,MPI::SUM,0);
    double acceptProb = exp(-netDeltaAction);
    bool acceptReject = RandomNumGenerator::getRand()>acceptProb;
    mpi->getWorkerComm().Bcast(&acceptReject,1,MPI::CHAR,0); 
    if (acceptReject) return false;
  } else {
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
    for (int imoving=0; imoving<nmoving; ++imoving) {
      Vec &bead(paths(movingIndex(imoving),islice));
      bead += displacement(imoving);
      bead = cell.pbc(bead);
    }
  }

  paths.setBuffers();

  return true;
}

AccRejEstimator* 
DisplaceMoveSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": moving " << nmoving
           << " " << particleChooser.getName();
  return accRejEst=new AccRejEstimator(longName.str(),1);
}
