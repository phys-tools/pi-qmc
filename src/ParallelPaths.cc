// $Id: ParallelPaths.cc,v 1.15 2007/07/28 02:42:07 jshumwa Exp $
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
#include "ParallelPaths.h"
#include "Beads.h"
#include "LinkSummable.h"
#include "Permutation.h"
#include "SuperCell.h"
#include <iostream>
#include <fstream>
#include "stats/MPIManager.h"
#include "BeadFactory.h"

ParallelPaths::ParallelPaths(int npart, int nslice, double tau,
    const SuperCell& cell, MPIManager &mpi, const BeadFactory &beadFactory)
  : Paths(npart,nslice,tau,cell), 
    iworker(mpi.getWorkerID()),
    nworker(mpi.getNWorker()),
    ifirst(nslice/nworker*iworker),
    nprocSlice(nslice/nworker+2+((iworker+1==nworker)?(nslice)%nworker:0)),
    beads(*beadFactory.getNewBeads(npart,nprocSlice)),
    buffer(*beadFactory.getNewBeads(npart,nprocSlice)),
    permutation(*new Permutation(npart)),
    inversePermutation(*new Permutation(npart)),
    mpi(mpi), npSlice(nworker) {
  if (mpi.isMain()) std::cout << "Creating parallel paths" << std::endl;
  std::cout <<  "workerID=" << mpi.getWorkerID() << std::endl;
  std::cout <<  "cloneID=" << mpi.getCloneID() << std::endl;
  for (int i=0; i<nworker; ++i)
    npSlice=nslice/nworker+1+((i+1==nworker)?nslice%nworker:0);
}

ParallelPaths::~ParallelPaths() {
  delete &beads; delete &buffer; 
  delete &permutation; delete &inversePermutation;
}

void ParallelPaths::sumOverLinks(LinkSummable& estimator) const {
  estimator.initCalc(nprocSlice-2,1+ifirst);
  for (int islice=1; islice<nprocSlice-1; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      estimator.handleLink(beads(ipart,islice-1), beads(ipart,islice),
                           ipart, (islice+ifirst)%nslice, *this);
    }
  }
  estimator.endCalc(nprocSlice-2);
}

Paths::Vec& 
ParallelPaths::operator()(const int ipart, const int islice) {
  return beads(ipart,(islice-ifirst+nslice)%nslice);
}

const Paths::Vec& 
ParallelPaths::operator()(const int ipart, const int islice) const {
  return beads(ipart,(islice-ifirst+nslice)%nslice);
}

Paths::Vec&
ParallelPaths::operator()(const int ipart, const int islice,
                                const int istep) {
  return beads(ipart,(islice-ifirst+istep+nslice)%nslice);
}

const void* ParallelPaths::getAuxBead(const int ipart, const int islice, 
    const int iaux) const {
  return beads.getAuxBead(ipart,(islice-ifirst+nslice)%nslice,iaux);
}

void* ParallelPaths::getAuxBead(const int ipart, const int islice, 
    const int iaux) {
  return beads.getAuxBead(ipart,(islice-ifirst+nslice)%nslice,iaux);
}

const Paths::Vec&
  ParallelPaths::operator()(const int ipart, const int islice,
                                  const int istep) const {
  return beads(ipart,(islice-ifirst+istep+nslice)%nslice);
}

Paths::Vec ParallelPaths::delta(const int ipart, const int islice,
                       const int istep) const {
  int jslice=(islice-ifirst+nslice)%nslice;
  Vec v = beads(ipart,jslice);
  v-=beads(ipart,jslice+istep);
  cell.pbc(v);
  return v;
}

void ParallelPaths::getBeads(const int ifirstSlice, 
                                   Beads<NDIM>& outBeads) const {
  int nsectionSlice=outBeads.getNSlice();
  int jfirstSlice=(ifirstSlice-ifirst+nslice)%nslice;
  for (int isectionSlice=0; isectionSlice<nsectionSlice; ++isectionSlice) {
    int islice=isectionSlice+jfirstSlice;
    for (int ipart=0; ipart<npart; ++ipart) {
      outBeads(ipart,isectionSlice)=beads(ipart,islice);
    }
  }
}

void ParallelPaths::getSlice(int islice, VArray& out) const {
  int jslice=(islice-ifirst+nslice)%nslice;
  for (int ipart=0; ipart<npart; ++ipart) out(ipart)=beads(ipart,jslice);
}

void ParallelPaths::putBeads(int ifirstSlice, const Beads<NDIM>& inBeads,
                           const Permutation& inPermutation) const {
  //First place the section beads back in the bead array.
  int jfirstSlice=(ifirstSlice-ifirst+nslice)%nslice;
  int nsectionSlice=inBeads.getNSlice();
  for (int isectionSlice=0; isectionSlice<nsectionSlice; ++isectionSlice) {
    int islice=isectionSlice+jfirstSlice;
    for (int ipart=0; ipart<npart; ++ipart) {
      beads(ipart,islice)=inBeads(ipart,isectionSlice);
    }
  }
  // Now permute the following beads.
  int jlastSlice=jfirstSlice+nsectionSlice;
  if (!inPermutation.isIdentity()){
    for (int islice=jlastSlice; islice<nprocSlice; ++islice) {
      // Swap paths.
      VArray buffer(npart);
      for (int i=0; i<npart; ++i) buffer(i)=beads(inPermutation[i],islice);
      for (int i=0; i<npart; ++i) beads(i,islice)=buffer(i);
    }
    permutation.prepend(inPermutation);
    inversePermutation.setToInverse(permutation);
  }
}

void ParallelPaths::putDoubleBeads(
          int ifirstSlice1,Beads<NDIM> &beads1, Permutation &p1,
          int ifirstSlice2,Beads<NDIM> &beads2, Permutation &p2) const {
  putBeads(ifirstSlice1,beads1,p1);
  putBeads(ifirstSlice2,beads2,p2);
}

void ParallelPaths::shift(const int ishift) {
#ifdef ENABLE_MPI
  // Transfer ishift beads between workers.
  int idest=(iworker+nworker-1)%nworker;
  int isrc=(iworker+1)%nworker;
  double *sendbuf=&(beads(0,0)[0]);
  double *recvbuf=&(buffer(0,nprocSlice-ishift-2)[0]);
  mpi.getWorkerComm().Sendrecv(
       sendbuf,(ishift+2)*NDIM*npart,MPI::DOUBLE,idest,1,
       recvbuf,(ishift+2)*NDIM*npart,MPI::DOUBLE,isrc,1);
  // Permute the ishift beads that we got from the next worker.
  for (int islice=nprocSlice-ishift-2; islice<nprocSlice; ++islice) {
    VArray vbuffer(npart);
    for (int i=0; i<npart; ++i) vbuffer(i)=buffer(permutation[i],islice);
    for (int i=0; i<npart; ++i) buffer(i,islice)=vbuffer(i);
  }
  // Shift the remaining beads that are in this worker.
  for (int j=0; j<nprocSlice-ishift-2; ++j) {
    for (int i=0; i<npart; ++i) {
      buffer(i,j)=beads(i,j+ishift);
    }
  }
  mpi.getWorkerComm().Barrier();
  cycleArrays(buffer.getCoordArray(),beads.getCoordArray());
#endif
} 

void ParallelPaths::setBuffers() {
  shift(0);
}
/*#ifdef ENABLE_MPI
  // Copy first two slices of beads to previous worker's buffers.
  int isrc=(iworker+1)%nworker;
  int idest=(iworker+nworker-1)%nworker;
  Vec *sendbuf=beads.getCoordArray().data();
  Vec *recvbuf=buffer.getCoordArray().data();
  mpi.getWorkerComm().Sendrecv(
       sendbuf,2*NDIM*npart,MPI::DOUBLE,idest,1,
       recvbuf,2*NDIM*npart,MPI::DOUBLE,isrc,1);
  // Permute buffered beads and put into last two slices of beads.
  for (int i=0; i<npart; ++i) {
    beads(i,nprocSlice-2)=buffer(permutation[i],0);
    beads(i,nprocSlice-1)=buffer(permutation[i],1);
  }
//#endif
} */

void ParallelPaths::clearPermutation() {
  permutation.reset(); inversePermutation.reset();
}
