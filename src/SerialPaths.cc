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
#include "SerialPaths.h"
#include "Beads.h"
#include "BeadFactory.h"
#include "LinkSummable.h"
#include "util/Permutation.h"
#include "util/SuperCell.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <blitz/tinyvec.h>

SerialPaths::SerialPaths(int npart, int nslice, double tau,
    const SuperCell& cell, const BeadFactory &beadFactory)
  : Paths(npart,nslice,tau,cell),
    beads(*beadFactory.getNewBeads(npart,nslice)),
    buffer1(*beadFactory.getNewBeads(npart,(nslice>1000)?1000:nslice)),
    buffer2(*beadFactory.getNewBeads(npart,(nslice>1000)?1000:nslice)),
    permutation(*new Permutation(npart)),
    inversePermutation(*new Permutation(npart)) {
  std::cout << "Creating serial paths with " 
            << nslice << " slices." << std::endl;
}

SerialPaths::~SerialPaths() {
  delete &beads; delete &buffer1; delete &buffer2;
  delete &permutation; delete &inversePermutation;
}

void SerialPaths::sumOverLinks(LinkSummable& estimator) const {
  estimator.initCalc(nslice,0);
  for (int ipart=0; ipart<npart; ++ipart) {
    estimator.handleLink(beads(inversePermutation[ipart],nslice-1),
                         beads(ipart,0), ipart, 0, *this);
  }
  for (int islice=1; islice<nslice; ++islice) {
    for (int ipart=0; ipart<npart; ++ipart) {
      estimator.handleLink(beads(ipart,islice-1), beads(ipart,islice),
                           ipart, islice, *this);
    }
  }
// Code to print out slices
//  for (int islice=0; islice<nslice; ++islice) {
//  std::cout << "islice="<<islice << ", " << beads(0,islice) << std::endl;
//  }
// end code to print out slices
  estimator.endCalc(nslice);
}

Paths::Vec& SerialPaths::operator()(const int ipart, const int islice) {
  return (islice>=0) ? beads(ipart,islice)
                     : beads(inversePermutation[ipart],nslice+islice);
}

const Paths::Vec& 
SerialPaths::operator()(const int ipart, const int islice) const {
  return (islice>=0) ? beads(ipart,islice)
                     : beads(inversePermutation[ipart],nslice+islice);
}

Paths::Vec&
SerialPaths::operator()(const int ipart, const int islice, const int istep) {
  int jslice=islice+istep;
  if (jslice<0) {
    return beads(inversePermutation[ipart],jslice+nslice);
  } else if (jslice>=nslice) {
    return beads(permutation[ipart],jslice-nslice);
  } else {
    return beads(ipart,jslice);
  }
}

const Paths::Vec&
  SerialPaths::operator()(const int ipart, const int islice, const int istep)
    const {
  int jslice=islice+istep;
  if (jslice<0) {
    return beads(inversePermutation[ipart],jslice+nslice);
  } else if (jslice>=nslice) {
    return beads(permutation[ipart],jslice-nslice);
  } else {
    return beads(ipart,jslice);
  }
}

Paths::Vec SerialPaths::delta(const int ipart, const int islice,
                       const int istep) const {
  Vec v = beads(ipart,islice);
  int jslice=islice+istep;
  if (jslice<0) {
    v-= beads(inversePermutation[ipart],jslice+nslice);
  } else if (jslice>=nslice) {
    v-= beads(permutation[ipart],jslice-nslice);
  } else {
    v-= beads(ipart,jslice);
  }
  cell.pbc(v);
  return v;
}

void SerialPaths::getBeads(const int ifirstSlice, Beads<NDIM>& outBeads) const {
  int nsectionSlice=outBeads.getNSlice();
  for (int isectionSlice=0; isectionSlice<nsectionSlice; ++isectionSlice) {
    int islice=isectionSlice+ifirstSlice;
    if (islice<0) {
      beads.copySlice(inversePermutation,nslice+islice,
                      outBeads,isectionSlice);
    } else if (islice<nslice) {
      beads.copySlice(islice,outBeads,isectionSlice);
    } else {
      beads.copySlice(permutation,islice-nslice,outBeads,isectionSlice);
    }
  }
}

void SerialPaths::getSlice(int islice, VArray& out) const {
  for (int ipart=0; ipart<npart; ++ipart) out(ipart)=beads(ipart,islice);
}

void SerialPaths::putBeads(int ifirstSlice, const Beads<NDIM>& inBeads,
                           const Permutation& inPermutation) const {
  SerialPaths::putBeads(ifirstSlice, inBeads, inPermutation,
                        0, inBeads.getNSlice());
}

void SerialPaths::putBeads(int ifirstSlice, const Beads<NDIM>& inBeads,
         const Permutation& inPermutation, const int ifirst,
         const int nmslice) const {
  //First place the section beads back in the bead array.
  int nsectionSlice=nmslice;
  for (int isectionSlice=0; isectionSlice<nsectionSlice; ++isectionSlice) {
    int islice=isectionSlice+ifirstSlice+ifirst;
    if (islice<0) {
      // Do nothing, since this slice should not have been moved. 
    } else if (islice<nslice) {
      inBeads.copySlice(isectionSlice+ifirst,beads,islice);
      //for (int ipart=0; ipart<npart; ++ipart) {
      //  beads(ipart,islice)=inBeads(ipart,isectionSlice);
      //}
    } else {
      islice%=nslice;
      inBeads.copySlice(isectionSlice+ifirst,beads,permutation,islice);
      //for (int ipart=0; ipart<npart; ++ipart) {
      //  beads(permutation[ipart],islice)=inBeads(ipart,isectionSlice);
      //}
    }
  }
  // Now permute the following beads.
  if (!inPermutation.isIdentity()) {
  int ifollowing=ifirstSlice+nsectionSlice;
  Permutation temp(inPermutation);
  if (ifollowing>nslice) {
    ifollowing-=nslice;
    temp.prepend(inversePermutation);
    temp.append(permutation);
  }
//BUG with aux beads!!!!
  for (int islice=ifollowing; islice<nslice; ++islice) {
    // Swap paths.
    VArray buffer(npart);
    for (int i=0; i<npart; ++i) buffer(i)=beads(temp[i],islice);
    for (int i=0; i<npart; ++i) beads(i,islice)=buffer(i);
  }
  permutation.prepend(temp);
  inversePermutation.setToInverse(permutation);
  }
}

void SerialPaths::putDoubleBeads(
          int ifirstSlice1,Beads<NDIM> &beads1, Permutation &p1,
          int ifirstSlice2,Beads<NDIM> &beads2, Permutation &p2) const {
  if (ifirstSlice1<nslice/2-1) {
    putBeads(ifirstSlice1,beads1,p1);
    Permutation p1inv(p1); p1inv.setToInverse(p1);
    p2.prepend(p1);
    p2.append(p1inv);
    beads2.permute(p1);
    putBeads(ifirstSlice2,beads2,p2);
  } else {
    putBeads(ifirstSlice2,beads2,p2);
    Permutation p2inv(p2); p2inv.setToInverse(p2);
    p1.prepend(p2);
    p1.append(p2inv);
    beads1.permute(p2);
    putBeads(ifirstSlice1,beads1,p1);
  }
}

const void* SerialPaths::getAuxBead(const int ipart, const int islice, 
    const int iaux) const {
  return beads.getAuxBead(ipart,islice,iaux);
}

void* SerialPaths::getAuxBead(const int ipart, const int islice, 
    const int iaux) {
  return beads.getAuxBead(ipart,islice,iaux);
}

void SerialPaths::shift(const int ishift) {
  if (ishift==0) return;
  Permutation p(npart);
  int nshift=(ishift>(buffer1.getNSlice())-1)?(buffer1.getNSlice()-1):ishift;
  int remainder=ishift-nshift;
  int islice=nslice-1;
  getBeads(islice,buffer1);
  for (islice-=nshift; islice>0; islice-=nshift) {
    getBeads(islice,buffer2); 
    putBeads(islice,buffer1,p,1,nshift);
    buffer1=buffer2;
  }
  putBeads(islice,buffer1,p,-islice,nshift+islice+1);
  // Recursively repeat if more shifting needed.
  if (remainder>0) shift(remainder);
}

void SerialPaths::clearPermutation() {
  permutation.reset(); inversePermutation.reset();
}
