 
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
#include "WritePaths.h"
#include "Paths.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "stats/MPIManager.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "SimulationInfo.h"

WritePaths::WritePaths(Paths& paths, const std::string& filename, const int dumpFreq, const int maxConfigs,  const bool writeMovie, const SimulationInfo& simInfo, MPIManager *mpi, const BeadFactory& beadFactory)  : filename(filename), paths(paths), mpi(mpi), beadFactory(beadFactory), dumpFreq(dumpFreq), maxConfigs(maxConfigs), writeMovie(writeMovie),simInfo(simInfo) {
  if (writeMovie) movieFile = new std::ofstream("pathMovie", std::ios::out);
}

void WritePaths::run() {
  
  //Write paths (movie) every dumpFreq (dumpMovieCounter) times
  static int dumpPathCounter=0;
  if (dumpPathCounter < dumpFreq) {
    dumpPathCounter++;
    return;
  }
  dumpPathCounter=0;
  static int dumpMovieCounter=0;
  dumpMovieCounter++;

  //Prepare file handlers 
  int workerID=(mpi)?mpi->getWorkerID():0;
  int nclone=(mpi)?mpi->getNClone():1;
  int cloneID=(mpi)?mpi->getCloneID():0;

   
  std::stringstream ext;
  if (nclone>1) ext << cloneID;
  std::ofstream *file=0;   
  //paths files
  Permutation perm(paths.getGlobalPermutation());   
  if (workerID==0){
    file = new std::ofstream((filename+ext.str()).c_str());
    *file <<"#Permutations ";
    for (int i=0; i<paths.getNPart(); i++) *file << perm[i]<< " ";
    *file << std::endl;
    *file << "#Path coordinates: " << paths.getNPart()
	  << " particles in " << NDIM << "-d" << std::endl;
  }
  //movie files
  if (cloneID==0 && workerID==0 && writeMovie) {
    if (dumpMovieCounter > maxConfigs) {
      (*movieFile).close();
      (*movieFile).open("pathMovie", std::ios::trunc);
      dumpMovieCounter=0;
    }
    *movieFile <<"#Permutations ";
    for (int i=0; i<paths.getNPart(); i++) *movieFile << perm[i]<< " ";
    *movieFile << std::endl;
    *movieFile << "# " << paths.getNPart()<<"   ";
      for (int ispec=0;ispec<simInfo.getNSpecies(); ispec++) *movieFile <<simInfo.getSpecies(ispec).count<<"  ";
      *movieFile <<std::endl;
  }

  /// Now determine how many chunked steps we need to write out all slices.
  int nslice=paths.getNSlice();
  int ifirst=paths.getLowestOwnedSlice(false);
  int imax=paths.getHighestOwnedSlice(false);
#ifdef ENABLE_MPI
  int nworker=(mpi)?mpi->getNWorker():1;
  if (mpi && nworker>1) {
    mpi->getWorkerComm().Bcast(&nslice,1,MPI::INT,0);
    mpi->getWorkerComm().Bcast(&ifirst,1,MPI::INT,0);
    mpi->getWorkerComm().Bcast(&imax,1,MPI::INT,0);
  }
#endif

  /*if (imax-ifirst==nslice) {
    // Use old method if no chunking.
    if (workerID==0){
      for (int ipart=0; ipart<paths.getNPart(); ++ipart) {
        Paths::Vec p = paths(ipart,paths.getNSlice()-1,1);
        for (int i=0; i<NDIM; ++i) *file << p[i] << " ";
      }
      *file << std::endl;
      for (int islice=paths.getNSlice()-1; islice>=0; --islice) {
        for (int ipart=0; ipart<paths.getNPart(); ++ipart) {
          Paths::Vec p = paths(ipart,islice);
          for (int i=0; i<NDIM; ++i) *file << p[i] << " ";
        }
        *file << std::endl;
      }
    }
  } else { */
    int nslicePerChunk = (imax-ifirst)/2;

    /// Loop over chunks.
    for (int ichunk=0; ichunk<ceil((nslice+1.)/nslicePerChunk); ++ichunk) {
      int lastSlice=(nslice+1-ifirst)-ichunk*nslicePerChunk;
      if (lastSlice>nslicePerChunk) lastSlice=nslicePerChunk;
      if (workerID==0){
        for (int islice=0; islice<lastSlice; ++islice) {
        //for (int islice=nslicePerChunk-1; 
        //         islice>=(nslicePerChunk-lastSlice)+ifirst; --islice) {
          for (int ipart=0; ipart<paths.getNPart(); ++ipart) {
            Paths::Vec p = paths(ipart,islice);
            for (int i=0; i<NDIM; ++i) *file << p[i] << " ";
	    if (writeMovie)  for (int i=0; i<NDIM; ++i) *movieFile << p[i] << " ";
          }
          *file << std::endl;
	  if (writeMovie) *movieFile << std::endl;
        }
      }
      if (ichunk+1<ceil((nslice+1.)/nslicePerChunk)) 
        paths.shift(nslicePerChunk);
    }
//  }

  if (workerID==0) delete file;
}

/*void WritePaths::run() {
  int workerID=(mpi)?mpi->getWorkerID():0;
  int nclone=(mpi)?mpi->getNClone():1;
  int nworker=(mpi)?mpi->getNWorker():1;
  int cloneID=(mpi)?mpi->getCloneID():0;
  std::ofstream *file=0;
  if (workerID==0){
    /// Add cloneID to name if there are clones.
    std::stringstream ext;
    if (nclone>1) ext << cloneID;
    std::cout << "Writting paths to file " << (filename+ext.str()) << std::endl;
    file = new std::ofstream((filename+ext.str()).c_str());
    *file << "#Path coordinates: " << paths.getNPart()
          << " particles in " << NDIM << "-d" << std::endl;
  }
  int npart=paths.getNPart();
  int nslice=paths.getNSlice();
  Beads<NDIM> &slice(*beadFactory.getNewBeads(npart,1));
  Permutation p(npart);
#ifdef ENABLE_MPI
  //if (mpi) mpi->getWorkerComm().Bcast(&slice(0,0),npart*NDIM,MPI::DOUBLE,0);
  if (paths.isProcessorSlice(0) && workerID!=0) {
     paths.getBeads(0,slice);
    if (mpi) mpi->getWorkerComm().Send(
                     &slice(0,0),npart*NDIM,MPI::DOUBLE,0,0);
  }
  if (workerID==0) {
    if (!paths.isProcessorSlice(0)) {
      if (mpi) mpi->getWorkerComm().Recv(
                       &slice(0,0),npart*NDIM,MPI::DOUBLE,MPI::ANY_SOURCE,0);
    } else {
       paths.getBeads(0,slice);
    }
  }
#endif
  if (workerID==0){
    for (int i=0; i<npart; ++i) {
      for (int idim=0; idim<NDIM; ++idim) *file << slice(i,0)[idim] << " ";
    }
    *file << std::endl;
  }
  for (int islice=0; islice<(nslice/bfactor)-1; ++islice) { 
    for (int ib=0; ib<bfactor; ++ib) { //Put the previous slice.
      int jslice=islice+ib*nfslice;
      if (paths.isProcessorSlice(jslice)) paths.putBeads(jslice,slice,p);
    }
    if (workerID==0){
    //if (!mpi || mpi->isMain()) {
      for (int i=0; i<npart; ++i) { //Read the next slice.
#if NDIM==3
        double x,y,z; *infile >> x >> y >> z;
        slice(i,0)=Beads<NDIM>::Vec(x,y,z);
#else
        double x,y; *infile >> x >> y;
        slice(i,0)=Beads<NDIM>::Vec(x,y);
#endif
      }
    }
#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(&slice(0,0),npart*NDIM,MPI::DOUBLE,0);
//  if (mpi) MPI::COMM_WORLD.Bcast(&slice(0,0),npart*NDIM,MPI::DOUBLE,0);
#endif
  }
  // Calculate  permuatation for last slice.
  if (workerID==0){
  //if (!mpi || mpi->isMain()) {
    Beads<NDIM> &wrapSlice(*beadFactory.getNewBeads(npart,1));
    for (int i=0; i<npart; ++i) {
#if NDIM==3
      double x,y,z; *infile >> x >> y >> z; 
      wrapSlice(i,0)=Beads<NDIM>::Vec(x,y,z);
#else
      double x,y; *infile >> x >> y; wrapSlice(i,0)=Beads<NDIM>::Vec(x,y);
#endif
    }
    for (int i=0; i<npart; ++i) {
      for (int j=0; j<npart; ++j) {
        if (sum(fabs(wrapSlice(i,0)-firstSlice(j,0)))<10e-6) {p[i]=j; }
      }
    }
    delete &wrapSlice;
  }   
#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(&p[0],npart,MPI::INT,0);
//  if (mpi) MPI::COMM_WORLD.Bcast(&p[0],npart,MPI::INT,0);
#endif
  for (int ib=bfactor-1; ib>=0; --ib) {
    int jslice=nfslice-1+ib*nfslice;
    if (paths.isProcessorSlice(jslice)) paths.putBeads(jslice,slice,p);
  }
  paths.setBuffers();
  delete file;
  if (workerID==0){
  //if (!mpi || mpi->isMain()) {
    std::cout << "Permuation read in is " << p << std::endl;
  } 
  delete &slice;
}
*/
 
