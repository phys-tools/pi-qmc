// $Id$
/*  Copyright (C) 2008 John B. Shumway, Jr.

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
#include "FreeEnergyEstimator.h"
#include "SimulationInfo.h"
#include "Paths.h"
#include "Permutation.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include <fstream>
#include "Species.h"
#include "EnumeratedModelState.h"
//sak
#include <iostream>
FreeEnergyEstimator::FreeEnergyEstimator(const SimulationInfo& simInfo,
  int nmodel, MPIManager *mpi)
  : BlitzArrayBlkdEst<1>("free_energy", "histogram/free-energy",
                         IVecN(nmodel), true), 
    mpi(mpi) {
  value = 0.;
  norm = 0;
}

FreeEnergyEstimator::~FreeEnergyEstimator() {
}

void FreeEnergyEstimator::evaluate(const Paths &paths) {
  const EnumeratedModelState *modelState 
    = dynamic_cast<const EnumeratedModelState*> (paths.getModelState());
  int imodel = modelState->getModelState();

  value(imodel) += 1.;
  norm += 1.;
}



void FreeEnergyEstimator::averageOverClones(const MPIManager* mpi) {
#ifdef ENABLE_MPI
  if (mpi && mpi->isCloneMain()) {
    int rank = mpi->getCloneComm().Get_rank();
    int size = mpi->getCloneComm().Get_size();
    if (size>1) {
      reset();
      if (rank==0) {
#if MPI_VERSION==2
        mpi->getCloneComm().Reduce(MPI::IN_PLACE,&norm,1,MPI::DOUBLE,
                                   MPI::SUM,0);
        mpi->getCloneComm().Reduce(MPI::IN_PLACE,value.data(),value.size(),
                                   MPI::FLOAT,MPI::SUM,0);
#else
        double nbuff;
        ArrayN vbuff(n);
        mpi->getCloneComm().Reduce(&norm,&nbuff,1,MPI::DOUBLE,
                                   MPI::SUM,0);
        mpi->getCloneComm().Reduce(value.data(),vbuff.data(),value.size(),
                                   MPI::FLOAT,MPI::SUM,0);
        norm=nbuff;
        value=vbuff;
#endif
      } else {
        mpi->getCloneComm().Reduce(&norm,NULL,1,MPI::DOUBLE,MPI::SUM,0);
        mpi->getCloneComm().Reduce(value.data(),NULL,value.size(),
                                   MPI::FLOAT,MPI::SUM,0);
      }
    }
  }
#endif
  // Next add value to accumvalue and accumvalue2.
  accumvalue += value/norm;
  if (hasErrorFlag) accumvalue2 += (value*value)/(norm*norm);
  accumnorm+=1.; 
  value=0.; norm=0;
  ++iblock; 
}
