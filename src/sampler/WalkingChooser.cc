//$Id$
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
#include "util/RandomNumGenerator.h"
#include "util/Permutation.h"
#include "MultiLevelSampler.h"
#include "Species.h"
#include "util/SuperCell.h"
#include "Beads.h"
#include "SimulationInfo.h"
#include "util/PeriodicGaussian.h"
#include <blitz/tinyvec-et.h>

WalkingChooser::WalkingChooser(const int nsize, const Species &species,
        const int nlevel, const SimulationInfo& simInfo)
:   PermutationChooser(nsize), SpeciesParticleChooser(species,nsize),
    t(npart,npart),
    cump(npart,npart),
    pg(NDIM),
    nsize(nsize),
    mass(species.mass),
    tau(simInfo.getTau()) {
    double alpha=mass/(2*tau*pow(2,nlevel));
    for (int idim=0; idim<NDIM; ++idim) {
        double length = (*simInfo.getSuperCell())[idim];
        pg(idim) = new PeriodicGaussian(alpha, length);
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
       if (index(ipart)==npart) index(ipart)=npart-1;
       //std :: cout <<"\n"<<"Begin : "<<index(0)<<"---->";
    } else { 
      // Select subsequent particle from t(k1,k2)/h(k1)
      double h=cump(index(ipart-1),npart-1);  
      double x=h*RandomNumGenerator::getRand();
      index(ipart)=iSearch(index(ipart-1),x);

      /*test JS and this method
	double xx=x;
	for (int jpart=0; jpart<npart-1; ++jpart) {
	xx-=t(index(ipart-1),jpart);
	if (xx<0) {
	if (index(ipart)==jpart) std :: cout <<"FOUND Agreement"<<index(ipart)<<" is equal"<<jpart<<std:: endl; 
	if (index(ipart)!=jpart) std :: cout <<"DIS DIS Agreement"<<index(ipart)<<" is NOT equal"<<jpart<<std:: endl; 
	break;}
	}
      */

      // Reject if repeated (not a nsize permutiaton cycle).
      for (int jpart=0; jpart<ipart; ++jpart) {
        if(index(ipart)==index(jpart)) return false;
      }
      // Accumulate inverse probabilities.
      revTranProb+=h/t(index(ipart-1),index(ipart-1));
      tranProb+=h/t(index(ipart-1),index(ipart));
      //std :: cout <<index(ipart)<<"---->";
    }
  }
  // Accumulate the inverse probability to connect the cycle.
  double h=cump(index(nsize-1),npart-1); 
  revTranProb+=h/t(index(nsize-1),index(nsize-1));
  tranProb+=h/t(index(nsize-1),index(0));
  double accept=revTranProb/tranProb;
  //std :: cout <<index(nsize-1)<<"END"<<std :: endl<<std :: endl;

  prob=1;
  for (int i=0; i<nsize; ++i) prob*=t(index(i),index(i))
                                   /t(index(i),index((i+1)%nsize));
  // Add ifirst to particle IDs to convert to absolute IDs.
  index += ifirst;
  return RandomNumGenerator::getRand()<accept;
}

/*
bool WalkingChooser::choosePermutation() {
  double tranProb=0, revTranProb=0;
  for (int ipart=0; ipart<nsize; ++ipart) {
    if (ipart== 0) {
       // Select first particle at random.
       index(ipart)=(int)(npart*RandomNumGenerator::getRand());
       if (index(ipart)==npart) index(ipart)=npart-1;
    } else { 
      // Select subsequent particle from t(k1,k2)/h(k1)
      double h=0;
      for (int jpart=0; jpart<npart; ++jpart) h+=t(index(ipart-1),jpart);
      double x=h*RandomNumGenerator::getRand();
      
      index(ipart)=npart-1;
      //int kpart; 
      for (int jpart=0; jpart<npart-1; ++jpart) {
      //int iorg=RandomNumGenerator::getRand()*(npart-1);
      //for (int jpart=iorg; jpart<(iorg+npart-1); ++jpart) {
      x-=t(index(ipart-1),jpart);
      if (x<0) {index(ipart)=jpart; break;}
      //kpart=jpart%(npart-1);//array[jpart];
      //x-=t(index(ipart-1),kpart);
      //if (x<0) {index(ipart)=kpart; break;}
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
  // Add ifirst to particle IDs to convert to absolute IDs. 
  index += ifirst;
  return RandomNumGenerator::getRand()<accept;
}
*/
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
              t(ipart,jpart) *= pg(idim)->evaluate(delta[idim]);
          }
      }
  }

 
  for (int ipart=0; ipart<npart; ++ipart) {
    cump(ipart,0)=t(ipart,0);
  }
  
  for (int ipart=0; ipart<npart; ++ipart) {
    for (int jpart=1; jpart<npart; ++jpart) {
      cump(ipart,jpart)=cump(ipart,jpart-1)+t(ipart,jpart);
      }
  }
}

int WalkingChooser::iSearch(int part, double x){
    
  if (x<cump(part,0))   return 0;
  
  if(cump(part,npart-2) < x) return npart-1;
  
  int ilo=0;
  int ihi=npart-2;
  int mid=0;
  for (int k=0;k<npart-1;k++){
    if(ihi==ilo+1) return ihi;
    mid=(ilo+ihi)/2;
    if(cump(part,mid)>x) {
      ihi=mid;
    } else{
      ilo=mid;
    }
  }
  return mid;
}

void WalkingChooser::setMLSampler(const MultiLevelSampler *mls) {
  multiLevelSampler=mls;
}
