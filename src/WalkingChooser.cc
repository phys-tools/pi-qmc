//$Id: WalkingChooser.cc,v 1.16 2006/10/18 17:08:19 jshumwa Exp $
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
#include "WalkingChooser.h"
#include "RandomNumGenerator.h"
#include "Permutation.h"
#include "MultiLevelSampler.h"
#include "Species.h"
#include "SuperCell.h"
#include "Beads.h"
#include "SimulationInfo.h"
#include "PeriodicGaussian.h"
#include <blitz/tinyvec-et.h>

WalkingChooser::WalkingChooser(const int nsize, const Species &species,
  const int nlevel, const SimulationInfo& simInfo)
  : PermutationChooser(nsize), SpeciesParticleChooser(species,nsize),
    nsize(nsize), t(npart,npart), mass(species.mass), 
    tau(simInfo.getTau()),
    pg(NDIM) {
  double alpha=mass/(2*tau*pow(2,nlevel));
  for (int idim=0; idim<NDIM; ++idim) {
    double l=(*simInfo.getSuperCell())[idim];
    int ngrid = (alpha*l*l>1) ? (int)(l*sqrt(alpha)*10) : 10;
    pg(idim)=new PeriodicGaussian(alpha,l,ngrid);
  }
  // Initialize permutation to an n-cycle.
  for (int i=0; i<nsize; ++i) (*permutation)[i]=(i+1)%nsize;
}

WalkingChooser::~WalkingChooser() {
  for (PGArray::iterator i=pg.begin();  i!=pg.end(); ++i) delete *i;
}

void WalkingChooser::chooseParticles() {}

bool WalkingChooser::choosePermutation() {
  double tranProb=0, revTranProb=0;
  for (int ipart=0; ipart<nsize; ++ipart) {
    if (ipart== 0) {
       // Select first particle at random.
       index(ipart)=(int)(npart*RandomNumGenerator::getRand());
       if (index(ipart)==npart) index(ipart)==npart-1;
    } else { 
      // Select subsequent particle from t(k1,k2)/h(k1)
      double h=0;
      for (int jpart=0; jpart<npart; ++jpart) h+=t(index(ipart-1),jpart);
      double x=h*RandomNumGenerator::getRand();
      index(ipart)=npart-1;
      for (int jpart=0; jpart<npart-1; ++jpart) {
        x-=t(index(ipart-1),jpart);
        if (x<0) {index(ipart)=jpart; break;}
      }
      // Reject if repeated (not a nsize permutiaton cycle).
      for (int jpart=0; jpart<ipart; ++jpart) {
        if(index(ipart)==index(jpart)) return false;
      }
      // Accumulate inverse probabilities.
      revTranProb+=h/t(index(ipart-1),index(ipart-1));
      tranProb+=h/t(index(ipart-1),index(ipart));
    }
  }
  // Accumulate the inverse probability to connect the cycle.
  double h=0;
  for (int jpart=0; jpart<npart; ++jpart) h+=t(index(nsize-1),jpart);
  revTranProb+=h/t(index(nsize-1),index(nsize-1));
  tranProb+=h/t(index(nsize-1),index(0));
  double accept=revTranProb/tranProb;
  prob=1;
  for (int i=0; i<nsize; ++i) prob*=t(index(i),index(i))
                                   /t(index(i),index((i+1)%nsize));
//std::cout << prob << std::endl;
  // Add ifirst to particle IDs to convert to absolute IDs.
  index += ifirst;
  return RandomNumGenerator::getRand()<accept;
}

void WalkingChooser::init() {
  // Setup the table of free particle propagotor values.
  const Beads<NDIM> &sectionBeads=multiLevelSampler->getSectionBeads();
  const SuperCell &cell=multiLevelSampler->getSuperCell();
  const int nslice=sectionBeads.getNSlice();
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int jpart=0; jpart<npart; ++jpart) {
      Vec delta=sectionBeads(ipart+ifirst,0);
      delta-=sectionBeads(jpart+ifirst,nslice-1);
      cell.pbc(delta);
      t(ipart,jpart)=1;
      for (int idim=0; idim<NDIM; ++idim) {
        t(ipart,jpart)*=(*pg(idim))(fabs(delta[idim]));
      }
    }
  }
}

void WalkingChooser::setMLSampler(const MultiLevelSampler *mls) {
  multiLevelSampler=mls;
}
