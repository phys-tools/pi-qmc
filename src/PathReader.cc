// $Id$
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
#include "PathReader.h"
#include "Paths.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Beads.h"
#include "BeadFactory.h"
#include "Permutation.h"
#include <blitz/tinyvec-et.h>
#include "stats/MPIManager.h"

void PathReader::run() {
  paths.clearPermutation();
  int workerID=(mpi)?mpi->getWorkerID():0;
  int nclone=(mpi)?mpi->getNClone():1;
  int cloneID=(mpi)?mpi->getCloneID():0;
  std::ifstream *infile=0;
  if (workerID==0){
    /// Add cloneID to name if there are clones.
    std::stringstream ext;
    if (nclone>1) ext << cloneID;
    std::cout << "Reading paths from file " << (filename+ext.str()) << std::endl;
    infile = new std::ifstream((filename+ext.str()).c_str());
  }
  int npart=paths.getNPart();
  int nslice=paths.getNSlice();
  int nfslice=nslice/bfactor;
  Beads<NDIM> &slice(*beadFactory.getNewBeads(npart,1));
  Permutation p(npart);   
  bool permutationsFlag = false;
  if (workerID==0){
    std::string temp; 
    *infile >> temp; 
    if (temp.compare("#Permutations") == 0) {
      std::cout << "Clone ID :: " << cloneID
         << " :: Found a Permutations array in restart file: " << std::endl;
      permutationsFlag = true;
      std :: cout << temp<<"  ";
      for (int i=0; i < npart; i++){
	*infile >>p[i];  
      }
      std :: cout<<p <<std :: endl;      
      getline(*infile,temp); 
      getline(*infile,temp); 
    }
    else{
      std::cout << "WARNING :: Did not find #Permutations in the restart file."
                << std::endl; 
      std::cout << "Automatically set up Permutations from paths."
                << std::endl; 
      getline(*infile,temp); 
    }
    for (int i=0; i<npart; ++i) {
      Beads<NDIM>::Vec v;
      for (int idim=0; idim<NDIM; ++idim) *infile >> v[idim];
      slice(i,0)=v;
    }
  }

#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(&slice(0,0),npart*NDIM,MPI::DOUBLE,0);
#endif

  Beads<NDIM> firstSlice(slice);  Permutation pidentity(npart); 
  for (int islice=0; islice<(nslice/bfactor)-1; ++islice) { 
    for (int ib=0; ib<bfactor; ++ib) { //Put the previous slice.
      int jslice=islice+ib*nfslice;
      if (paths.isOwnedSlice(jslice)) paths.putBeads(jslice,slice,pidentity);
    }
    if (workerID==0){
    //if (!mpi || mpi->isMain()) {
      for (int i=0; i<npart; ++i) { //Read the next slice.
        Beads<NDIM>::Vec v;
        for (int idim=0; idim<NDIM; ++idim) *infile >> v[idim];
        slice(i,0)=v;
      }
    }
#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(&slice(0,0),npart*NDIM,MPI::DOUBLE,0);
#endif
  }

  
  if (!permutationsFlag) {
    // Calculate  permuatation for last slice.
    if (workerID==0){
      Beads<NDIM> &wrapSlice(*beadFactory.getNewBeads(npart,1));
      for (int i=0; i<npart; ++i) {
	Beads<NDIM>::Vec v;
	for (int idim=0; idim<NDIM; ++idim) *infile >> v[idim];
	wrapSlice(i,0)=v;
      }
      for (int i=0; i<npart; ++i) {
	for (int j=0; j<npart; ++j) {
	  if (sum(fabs(wrapSlice(i,0)-firstSlice(j,0)))<10e-6) {p[i]=j; }
	}
      }
      delete &wrapSlice;
    }   
  }


#ifdef ENABLE_MPI
    if (mpi) mpi->getWorkerComm().Bcast(&p[0],npart,MPI::INT,0);
#endif
  for (int ib=bfactor-1; ib>=0; --ib) {
    int jslice=nfslice-1+ib*nfslice;
    if (paths.isOwnedSlice(jslice)) paths.putBeads(jslice,slice,p);
  }
  paths.setBuffers();
  delete infile;
  if (workerID==0 && !permutationsFlag){
    std::cout << "Clone ID :: " << cloneID << " :: Permuation read in is "
              << p << std::endl; 
  }
  delete &slice;


  /*compare p amd perm
    for (int i=0; i<npart; ++i) {
      if(p[i] != perm[i] && permutationsFlag) std :: cout <<"ERROR :: clone "<<cloneID<<",  worker "<<workerID<<", discrepancy found for permutation pickup at part "<<i <<": p[i]="<<p[i]<<"perm[i]="<<perm[i]<<std :: endl<<std::flush;
    }
*/
}
 
