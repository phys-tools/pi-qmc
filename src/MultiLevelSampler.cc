// $Id$
/*  Copyright (C) 2004-2006,2009 John B. Shumway, Jr.

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
#include "MultiLevelSampler.h"
#include "stats/AccRejEstimator.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "Action.h"
#include "Mover.h"
#include "RandomNumGenerator.h"
#include "Paths.h"
#include "Permutation.h"
#include "PermutationChooser.h"
#include "ParticleChooser.h"
#include "SimpleParticleChooser.h"
#include "RandomPermutationChooser.h"
#include "SectionChooser.h"
#include <sstream>
#include <string>

MultiLevelSampler::MultiLevelSampler(int nmoving, Paths& paths,
  SectionChooser &sectionChooser, ParticleChooser& particleChooser,
  PermutationChooser& permutationChooser, Mover& mover, Action* action,
  const int nrepeat, const BeadFactory& beadFactory)
  : nlevel(sectionChooser.getNLevel()),
    nmoving(nmoving), sectionBeads(&sectionChooser.getBeads()),
    sectionPermutation(&sectionChooser.getPermutation()),
    movingBeads(beadFactory.getNewBeads(nmoving,sectionBeads->getNSlice())),
    mover(mover), cell(paths.getSuperCell()), action(action),
    movingIndex(new IArray(nmoving)),
    identityIndex(nmoving), pMovingIndex(nmoving),
    particleChooser(particleChooser),
    permutationChooser(permutationChooser),
    sectionChooser(sectionChooser), paths(paths), accRejEst(0),
    nrepeat(nrepeat) {
  for (int i=0; i<nmoving; ++i) (*movingIndex)(i)=identityIndex(i)=i;
}

MultiLevelSampler::~MultiLevelSampler() {
  delete movingBeads;
  delete movingIndex;
  delete &particleChooser;
  delete &permutationChooser;
}

void MultiLevelSampler::run() {
  // Select particles to move and the permuation.
  permutationChooser.init();
  for (int irepeat=0; irepeat<nrepeat; ++irepeat) {
    bool isNewPerm=permutationChooser.choosePermutation();
    if (isNewPerm) {
      Permutation permutation(permutationChooser.getPermutation());
      particleChooser.chooseParticles();
      double lnTranProb=permutationChooser.getLnTranProb();
      for (int i=0; i<nmoving; ++i) (*movingIndex)(i)=particleChooser[i];
      // Copy old coordinate endpoint to the moving coordinate endpoints.
      const int  nsectionSlice=movingBeads->getNSlice();
      for (int imoving=0; imoving<nmoving; ++imoving) {
        pMovingIndex(imoving)=(*movingIndex)(permutation[imoving]);
      }

     

      sectionBeads->copySlice(*movingIndex,0,*movingBeads,identityIndex,0);
       sectionBeads->copySlice(pMovingIndex,nsectionSlice-1,
                 *movingBeads,identityIndex,nsectionSlice-1);
//        (*movingBeads)(imoving,0)=(*sectionBeads)(particleChooser[imoving],0);
//        (*movingBeads)(imoving,nsectionSlice-1)
//          =(*sectionBeads)(particleChooser[permutation[imoving] ],
//                        nsectionSlice-1);
//      }
      if (tryMove(lnTranProb) && irepeat<nrepeat-1)  permutationChooser.init();
    }
  }
}

bool MultiLevelSampler::tryMove(double initialLnTranProb) {
  double oldDeltaAction=initialLnTranProb;
  for (int ilevel=nlevel; ilevel>=0; --ilevel) {
    if (accRejEst) accRejEst->tryingMove(ilevel);
    // Do the trial move for this level.
    double lnTranProb=mover.makeMove(*this,ilevel);
    // Evaluate the change in action for this level.
    double deltaAction=(action==0)?0:action->getActionDifference(*this,ilevel);
    double acceptProb=exp(lnTranProb-deltaAction+oldDeltaAction);
    //std::cout << acceptProb << ", " << lnTranProb << ",  " << ilevel << ", "
    //          << -deltaAction  << ", " << oldDeltaAction << std::endl;
    if (RandomNumGenerator::getRand()>acceptProb) return false;
    oldDeltaAction=deltaAction;
    if (accRejEst) accRejEst->moveAccepted(ilevel);
  }
  // Move accepted.
  action->acceptLastMove();
  // Put moved beads in section beads.
  for (int islice=0; islice<sectionBeads->getNSlice(); ++islice) {
    movingBeads->copySlice(identityIndex,islice,
              *sectionBeads,*movingIndex,islice);
//    for (int i=0; i<nmoving; ++i) {
//      (*sectionBeads)(particleChooser[i],islice)=(*movingBeads)(i,islice);
//    }
  }
  // Append the current permutation to section permutation.
  const Permutation& perm(permutationChooser.getPermutation());
  Permutation temp(perm);
  for (int i=0; i<nmoving; ++i) {
    temp[i]= (*sectionPermutation)[particleChooser[i]];
  }
  for (int i=0; i<nmoving; ++i) {
    (*sectionPermutation)[particleChooser[i]]=temp[perm[i]];
  }
  return true;
}

void MultiLevelSampler::setAction(Action* act, const int level) {action=act;}

AccRejEstimator* 
MultiLevelSampler::getAccRejEstimator(const std::string& name) {
  std::ostringstream longName;
  longName << name << ": level " << nlevel << ", moving " << nmoving
           << " " << particleChooser.getName();
  return accRejEst=new AccRejEstimator(longName.str(),nlevel+1);
}
