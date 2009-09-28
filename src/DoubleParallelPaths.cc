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
#include "DoubleParallelPaths.h"
#include "Beads.h"
#include "LinkSummable.h"
#include "Permutation.h"
#include "SuperCell.h"
#include <iostream>
#include <fstream>
#include "stats/MPIManager.h"
#include "BeadFactory.h"

DoubleParallelPaths::DoubleParallelPaths(int npart, int nslice, double tau,
    const SuperCell& cell, MPIManager &mpi, const BeadFactory &beadFactory)
  : Paths(npart,nslice,tau,cell), 
    iworker(mpi.getWorkerID()),
    nworker(mpi.getNWorker()),
    ifirst(nslice/2/nworker*iworker-1),
    nprocSlice(nslice/2/nworker+2+((iworker+1==nworker)?(nslice/2)%nworker:0)),
    beads1(*beadFactory.getNewBeads(npart,nprocSlice)),
    beads2(*beadFactory.getNewBeads(npart,nprocSlice)),
    buffer1(*beadFactory.getNewBeads(npart,nprocSlice)),
    buffer2(*beadFactory.getNewBeads(npart,nprocSlice)),
    permutation1(*new Permutation(npart)), globalPermutation(*new Permutation(npart)),
    permutation2(*new Permutation(npart)),  
    inversePermutation1(*new Permutation(npart)),
    inversePermutation2(*new Permutation(npart)),
    mpi(mpi), npSlice(nworker) {
  if (mpi.isMain()) std::cout << "Creating double parallel paths" << std::endl;
  for (int i=0; i<nworker; ++i)
    npSlice=nslice/2/nworker+1+((i+1==nworker)?(nslice/2)%nworker:0);
  std :: cout <<"CloneID :: "<< mpi.getCloneID()  <<" :: Creating DoubleParallelPaths on  WorkerID "<< mpi.getWorkerID()<<" with nprocSlice "<<nprocSlice<<std ::endl;
}

DoubleParallelPaths::~DoubleParallelPaths() {
  delete &beads1; delete &beads2; delete &buffer1; delete &buffer2;
  delete &permutation1; delete &inversePermutation1;
  delete &permutation2; delete &inversePermutation2;
}

void DoubleParallelPaths::sumOverLinks(LinkSummable& estimator) const {
  estimator.initCalc(2*(nprocSlice-2),1+ifirst);
  for (int islice=1; islice<nprocSlice-1; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      estimator.handleLink(beads1(ipart,islice-1), beads1(ipart,islice),
                           ipart, (islice+ifirst)%nslice, *this);
    }
  }
  for (int islice=1; islice<nprocSlice-1; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      estimator.handleLink(beads2(ipart,islice-1), beads2(ipart,islice),
                           ipart, (islice+ifirst+nslice/2)%nslice, *this);
    }
  }
  estimator.endCalc(2*(nprocSlice-2));
}

Paths::Vec& 
DoubleParallelPaths::operator()(const int ipart, const int islice) {
  return ((islice-ifirst+nslice)%nslice<nprocSlice)
           ?beads1(ipart,(islice-ifirst+nslice)%nslice)
           :beads2(ipart,islice-ifirst-nslice/2);
}

const Paths::Vec& 
DoubleParallelPaths::operator()(const int ipart, const int islice) const {
  return ((islice-ifirst+nslice)%nslice<nprocSlice)
           ?beads1(ipart,(islice-ifirst+nslice)%nslice)
           :beads2(ipart,islice-ifirst-nslice/2);
}

Paths::Vec&
DoubleParallelPaths::operator()(const int ipart, const int islice,
                                const int istep) {
  return ((islice-ifirst+nslice)%nslice<nprocSlice)
           ?beads1(ipart,(islice-ifirst+istep+nslice)%nslice)
           :beads2(ipart,(islice-ifirst+istep+nslice/2)%nslice);
}

const void* DoubleParallelPaths::getAuxBead(const int ipart, const int islice, 
    const int iaux) const {
  return ((islice-ifirst+nslice)%nslice<nprocSlice)
         ?beads1.getAuxBead(ipart,(islice-ifirst+nslice)%nslice,iaux)
         :beads2.getAuxBead(ipart,(islice-ifirst+nslice/2)%nslice,iaux);
}

void* DoubleParallelPaths::getAuxBead(const int ipart, const int islice, 
    const int iaux) {
  return ((islice-ifirst+nslice)%nslice<nprocSlice)
         ?beads1.getAuxBead(ipart,(islice-ifirst+nslice)%nslice,iaux)
         :beads2.getAuxBead(ipart,(islice-ifirst+nslice/2)%nslice,iaux);
}

const Paths::Vec&
  DoubleParallelPaths::operator()(const int ipart, const int islice,
                                  const int istep) const {
  return ((islice-ifirst+nslice)%nslice<nprocSlice)
           ?beads1(ipart,(islice-ifirst+istep+nslice)%nslice)
           :beads2(ipart,(islice-ifirst+istep+nslice/2)%nslice);
}

Paths::Vec DoubleParallelPaths::delta(const int ipart, const int islice,
                       const int istep) const {
  Beads<NDIM>& beads(((islice-ifirst)%nslice<nprocSlice)?beads1:beads2);
  int jslice=(islice-ifirst)%(nslice/2);
  Vec v = beads(ipart,jslice);
  v-=beads(ipart,jslice+istep);
  cell.pbc(v);
  return v;
}

void DoubleParallelPaths::getBeads(const int ifirstSlice, 
                                   Beads<NDIM>& outBeads) const {
  int nsectionSlice=outBeads.getNSlice();
  Beads<NDIM>& beads(((ifirstSlice-ifirst+nslice)%nslice<nprocSlice)?beads1:beads2);
  int jfirstSlice=(ifirstSlice-ifirst)%(nslice/2);
  for (int isectionSlice=0; isectionSlice<nsectionSlice; ++isectionSlice) {
    int islice=isectionSlice+jfirstSlice;
    for (int ipart=0; ipart<npart; ++ipart) {
      outBeads(ipart,isectionSlice)=beads(ipart,islice);
    }
  }
}

void DoubleParallelPaths::getSlice(int islice, VArray& out) const {
  Beads<NDIM>& beads(((islice-ifirst+nslice)%nslice<nprocSlice)?beads1:beads2);
  int jslice=(islice-ifirst)%(nslice/2);
  for (int ipart=0; ipart<npart; ++ipart) out(ipart)=beads(ipart,jslice);
}

void DoubleParallelPaths::putBeads(int ifirstSlice, const Beads<NDIM>& inBeads,
                           const Permutation& inPermutation) const {
  //First place the section beads back in the bead array.
  bool first = (((ifirstSlice-ifirst+nslice)%nslice)<nprocSlice);
  Beads<NDIM>& beads(first?beads1:beads2);
  Permutation& permutation(first?permutation1:permutation2);
  Permutation& inversePermutation(first
                                  ?inversePermutation1:inversePermutation2);
  int jfirstSlice=(ifirstSlice-ifirst)%(nslice/2);
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

void DoubleParallelPaths::putDoubleBeads(
          int ifirstSlice1,Beads<NDIM> &beads1, Permutation &p1,
          int ifirstSlice2,Beads<NDIM> &beads2, Permutation &p2) const {
  putBeads(ifirstSlice1,beads1,p1);
  putBeads(ifirstSlice2,beads2,p2);
}

const Permutation& DoubleParallelPaths::getGlobalPermutation() const {
  globalPermutation = permutation1;
  //return permutation1;
#ifdef ENABLE_MPI
  Permutation recvp1(npart); 
  Permutation recvp2(npart);
  Permutation perm2(permutation2);
 
  // get globalPerm
  if (iworker ==0) {
    for (int src =1; src < nworker; ++src){
  
      mpi.getWorkerComm().Recv(&recvp1[0], npart,MPI::INT,src,1);
      globalPermutation.append(recvp1);
    }
  } else {
    mpi.getWorkerComm().Send(&(permutation1[0]), npart,MPI::INT,0,1);
  }

  //get perm2
  if (iworker ==0) {
    for (int src =1; src < nworker; ++src){
      
      mpi.getWorkerComm().Recv(&recvp2[0], npart,MPI::INT,src,2);
      perm2.append(recvp2);
    }
  } else {
    mpi.getWorkerComm().Send(&(permutation2[0]), npart,MPI::INT,0,2);
  }   
  
  if (iworker ==0) globalPermutation.append(perm2);
  mpi.getWorkerComm().Bcast(&globalPermutation[0],npart,MPI::INT,0);
  
#endif   
  return globalPermutation;
}


void DoubleParallelPaths::shift(const int ishift) {
#ifdef ENABLE_MPI
  // Transfer ishift beads between workers.
  int idest=(iworker+nworker-1)%nworker;
  int isrc=(iworker+1)%nworker;
  double *sendbuf1=&(beads1(0,0)[0]);
  double *recvbuf1=(iworker<nworker-1)?&(buffer1(0,nprocSlice-ishift-2)[0])
                                      :&(buffer2(0,nprocSlice-ishift-2)[0]);
  mpi.getWorkerComm().Sendrecv(
       sendbuf1,(ishift+2)*NDIM*npart,MPI::DOUBLE,idest,1,
       recvbuf1,(ishift+2)*NDIM*npart,MPI::DOUBLE,isrc,1);
  double *sendbuf2=&(beads2(0,0)[0]);
  double *recvbuf2=(iworker<nworker-1)?&(buffer2(0,nprocSlice-ishift-2)[0])
                                      :&(buffer1(0,nprocSlice-ishift-2)[0]);
  mpi.getWorkerComm().Sendrecv(
       sendbuf2,(ishift+2)*NDIM*npart,MPI::DOUBLE,idest,2,
       recvbuf2,(ishift+2)*NDIM*npart,MPI::DOUBLE,isrc,2);
  // Permute the ishift beads that we got from the next worker.
  for (int islice=nprocSlice-ishift-2; islice<nprocSlice; ++islice) {
    VArray buffer(npart);
    for (int i=0; i<npart; ++i) buffer(i)=buffer1(permutation1[i],islice);
    for (int i=0; i<npart; ++i) buffer1(i,islice)=buffer(i);
    for (int i=0; i<npart; ++i) buffer(i)=buffer2(permutation2[i],islice);
    for (int i=0; i<npart; ++i) buffer2(i,islice)=buffer(i);
  }
  // Shift the remaining beads that are in this worker.
  for (int j=0; j<nprocSlice-ishift-2; ++j) {
    for (int i=0; i<npart; ++i) {
      buffer1(i,j)=beads1(i,j+ishift);
      buffer2(i,j)=beads2(i,j+ishift);
    }
  }
  mpi.getWorkerComm().Barrier();
  cycleArrays(buffer1.getCoordArray(),beads1.getCoordArray());
  cycleArrays(buffer2.getCoordArray(),beads2.getCoordArray());
#endif
} 

void DoubleParallelPaths::setBuffers() {
#ifdef ENABLE_MPI
  Beads<NDIM> buffer1(npart,2), buffer2(npart,2);
  // Copy first two slices of beads to previous worker's buffers.
  int isrc=(iworker+1)%nworker;
  int idest=(iworker+nworker-1)%nworker;
  Vec *sendbuf=beads1.getCoordArray().data();
  Vec *recvbuf=((iworker<nworker-1)?buffer1:buffer2).getCoordArray().data();
  mpi.getWorkerComm().Sendrecv(
       sendbuf,2*NDIM*npart,MPI::DOUBLE,idest,1,
       recvbuf,2*NDIM*npart,MPI::DOUBLE,isrc,1);
  sendbuf=beads2.getCoordArray().data();
  recvbuf=((iworker<nworker-1)?buffer2:buffer1).getCoordArray().data();
  mpi.getWorkerComm().Sendrecv(
       sendbuf,2*NDIM*npart,MPI::DOUBLE,idest,2,
       recvbuf,2*NDIM*npart,MPI::DOUBLE,isrc,2);
  // Permute buffered beads and put into last two slices of beads.
  for (int i=0; i<npart; ++i) {
    beads1(i,nprocSlice-2)=buffer1(permutation1[i],0);
    beads1(i,nprocSlice-1)=buffer1(permutation1[i],1);
    beads2(i,nprocSlice-2)=buffer2(permutation2[i],0);
    beads2(i,nprocSlice-1)=buffer2(permutation2[i],1);
  }
#endif
}

void DoubleParallelPaths::clearPermutation() {
  permutation1.reset(); inversePermutation1.reset();
  permutation2.reset(); inversePermutation2.reset();
}
