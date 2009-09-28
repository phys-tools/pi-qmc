// $Id: UniformMover.cc 22 2009-03-06 13:52:07Z john.shumwayjr $
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
#include "UniformMover.h"
#include "Beads.h"
#include "DisplaceMoveSampler.h"
#include "DoubleDisplaceMoveSampler.h"
#include "RandomNumGenerator.h"
//#include <blitz/tinyvec.h>
#include "SuperCell.h"
#include "SimulationInfo.h"
#include <cmath>


#include <sstream>
#include <string>
UniformMover::UniformMover(const MPIManager* mpi)
  : mpi(mpi) {
}

UniformMover::~UniformMover() {
 }

double UniformMover::makeMove(DisplaceMoveSampler& sampler) {
  // typedef blitz::TinyVector<double,NDIM> Vec;
  Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const double dist = sampler.getDist();
  const int nSlice = sampler.getNSlice();
  const int ifirstSlice = 0; 
   IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  Array uniRand(nMoving*NDIM);
  RandomNumGenerator::makeRand(uniRand); // fix::makerand should take a vec. also.  
#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    mpi->getWorkerComm().Bcast(&uniRand(0),nMoving*NDIM,MPI::DOUBLE,0);
  }
#endif
  
  // Calculate the new position.
  double * dr = new double[NDIM]; 
  for (int i=0; i<NDIM; i++){
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      dr[i] =  dist*(uniRand(NDIM*iMoving+i)-0.5);	
      for (int islice=0; islice<nSlice; islice++) {
	movingBeads(iMoving,islice)[i] += dr[i];
	cell.pbc(movingBeads(iMoving,islice));
      }
    }
  }

  delete [] dr;
  return 0; 
}


double UniformMover::makeMove(DoubleDisplaceMoveSampler& sampler) {
  // typedef blitz::TinyVector<double,NDIM> Vec;
  // Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  Beads<NDIM>& movingBeads1=sampler.getMovingBeads(1);
  Beads<NDIM>& movingBeads2=sampler.getMovingBeads(2); 
  const SuperCell& cell=sampler.getSuperCell();
  const double dist = sampler.getDist();
  const int nSlice = sampler.getNSlice(); 
 
  IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  Array uniRand(nMoving*NDIM);// fix::makerand should take a vec. also.  
  RandomNumGenerator::makeRand(uniRand); 
#ifdef ENABLE_MPI
  if (mpi && (mpi->getNWorker())>1) {
    mpi->getWorkerComm().Bcast(&uniRand(0),nMoving*NDIM,MPI::DOUBLE,0);
  }
#endif
 
  double * dr = new double[NDIM]; 
  for (int i=0; i<NDIM; i++){ 
    for (int iMoving=0; iMoving<nMoving; ++iMoving) {
      dr[i] =  dist*(uniRand(NDIM*iMoving+i)-0.5);	
      for (int islice=0; islice<nSlice; islice++) {
	//movingBeads(iMoving,islice)[i] += dr[i];
	movingBeads1(iMoving,islice)[i] += dr[i];
	cell.pbc(movingBeads1(iMoving,islice));
	
	movingBeads2(iMoving,islice)[i] += dr[i];
	cell.pbc(movingBeads2(iMoving,islice));
      }
    }
  }
  
  delete [] dr;
  return 0; 
}
