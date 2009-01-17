// $Id: CubicLattice.cc,v 1.16 2007/11/26 00:00:14 jshumwa Exp $
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
#include "RandomNumGenerator.h"
#include "stats/MPIManager.h"
#include "CubicLattice.h"
#include "Paths.h"
#include "Species.h"
#include <blitz/tinyvec-et.h>
#include "SuperCell.h"

CubicLattice::CubicLattice(Paths& paths, double a, double scatter,
     const IVec nmax, const Vec r0, MPIManager *mpi)
   : paths(paths), a(a), scatter(scatter), nmax(nmax), r0(r0), ifirst(0),
     npart(paths.getNPart()), mpi(mpi) {
}

CubicLattice::CubicLattice(Paths& paths, double a, double scatter,
     const IVec nmax, const Vec r0, const Species& species, MPIManager *mpi) 
   : paths(paths), a(a), scatter(scatter), nmax(nmax), r0(r0), 
     ifirst(species.ifirst), npart(species.count), mpi(mpi) {
}
		   

void CubicLattice::run() {
  int nslice=paths.getNSlice();
  int nsites=product(nmax);
  // Determine number of particles per site. Place equal numbers
  // of particles on each site, then distribute the remainder randomly.
  IArray n(nsites); 
  n=npart/nsites; 
  if (!mpi || mpi->isCloneMain()) {
    for (int ipart=0; ipart<npart%nsites;)  {
      int i=(int)(nsites*RandomNumGenerator::getRand());
      if (n(i)==npart/nsites) {++n(i); ++ipart;}
    }
  }
#ifdef ENABLE_MPI
  if (mpi) mpi->getWorkerComm().Bcast(n.data(),nsites,MPI::INT,0);
#endif
  // Now place the particles on the lattice.
  int ifirstSlice=paths.getLowestSampleSlice(0,false);
  int ilastSlice=paths.getHighestStoredSlice(0,false);
  SuperCell cell=paths.getSuperCell();
  if (!mpi || mpi->isCloneMain()) {
  int ipart=ifirst; int isite=0;
  for (int i=0; i<nmax[0]; ++i) {
#if NDIM==2
    for (int j=0; j<nmax[1]; ++j) {
#endif
#if NDIM==3
    for (int j=0; j<nmax[1]; ++j) {
      for (int k=0; k<nmax[2]; ++k) {
#endif
#if NDIM==4
    for (int j=0; j<nmax[1]; ++j) {
      for (int k=0; k<nmax[2]; ++k) {
       for (int l=0; l<nmax[3]; ++l) {
#endif
        for (int ii=0; ii<n(isite); ++ii) {
          paths(ipart,ifirstSlice)
#if NDIM==4
         =Vec(a*(i-(nmax[0]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)),
              a*(j-(nmax[1]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)),
              a*(k-(nmax[2]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)),
              a*(l-(nmax[3]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)) )
#endif
#if NDIM==3
         =Vec(a*(i-(nmax[0]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)),
              a*(j-(nmax[1]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)),
              a*(k-(nmax[2]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)))
#endif
#if NDIM==2
         =Vec(a*(i-(nmax[0]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)),
              a*(j-(nmax[1]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)))
#endif
#if NDIM==1
         =Vec(a*(i-(nmax[0]-1)/2.+scatter*(RandomNumGenerator::getRand()-0.5)))
#endif
              +r0;
	  cell.pbc(paths(ipart,ifirstSlice));
	  ++ipart;  
        }
        ++isite;
#if NDIM==2
       }
#endif	
#if NDIM==3
      }
    }
#endif	
#if NDIM==4
       }
      } 
    }
#endif
  }
  }
  // Copy first slice to workers.
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Bcast(&paths(ifirst,ifirstSlice), npart*NDIM,
                               MPI::DOUBLE,0);
  } 
#endif
  // Copy to other slices.
  for (int islice=ifirstSlice+1; islice<ilastSlice; ++islice) {
    for (int ipart=ifirst; ipart<ifirst+npart; ++ipart) {
      paths(ipart,islice)=paths(ipart,ifirstSlice);
    } 
  }
  if (paths.isDouble()) {
  for (int islice=ifirstSlice; islice<ilastSlice; ++islice) {
    for (int ipart=ifirst; ipart<ifirst+npart; ++ipart) {
      paths(ipart,islice+nslice/2)=paths(ipart,islice);
    } 
  }
  }
}
