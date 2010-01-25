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
#include "PermutationEstimator.h"
#include "SimulationInfo.h"
#include "Paths.h"
#include "Permutation.h"
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"
#include <fstream>
#include "Species.h"
//sak
#include <iostream>
PermutationEstimator::PermutationEstimator(const SimulationInfo& simInfo, const std::string& name,
		       const Species &s1, MPIManager *mpi)
  : BlitzArrayBlkdEst<1>(name, IVecN(s1.count), true), 
    ifirst(s1.ifirst), nipart(s1.count), npart(simInfo.getNPart()),
    mpi(mpi) {
  value = 0.;
  norm = 0;
  visited = new bool[nipart];
  for (int i=0; i<nipart; ++i)   visited[i] = false;
}

PermutationEstimator::~PermutationEstimator() {
  delete [] visited;
}

void PermutationEstimator::evaluate(const Paths &paths) {
  const Permutation &perm(paths.getGlobalPermutation());
  //std :: cout << "SAK :: IN PERMUTATION ESTIMATOR"<<std :: endl;
  //std :: cout << perm<<std :: endl;

  for (int k=0; k<nipart; k++){
    int i = ifirst + k;
    if (!visited[k]){
      int cnt = 1;
      visited[k] = true;
      int j = i;
      while(perm[j] !=i){
	cnt++;	
	j = perm[j];
	if ((j-ifirst) > (ifirst+nipart)){
	  std :: cout<< "*** WARNING in permutationEstimator"<<std ::endl;
	  break;
	}
	visited[j-ifirst]= true;
      }
      value(cnt) += 1.;
    }
  }
  for (int i=0; i<nipart; ++i)   visited[i] = false;
  norm += 1.;
  
}



void PermutationEstimator::averageOverClones(const MPIManager* mpi) {
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
