// $Id$
/*  Copyright (C) 2004-2009 John B. Shumway, Jr.

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
#include "EwaldImagePairAction.h"
#include "MultiLevelSampler.h"
#include "Beads.h"
#include "Paths.h"
#include "SuperCell.h"
#include "SimulationInfo.h"
#include <vector>
#include <fstream>
#include <blitz/tinyvec-et.h>
#include "DisplaceMoveSampler.h"

EwaldImagePairAction::EwaldImagePairAction(const Species& s1, const Species& s2,
  const EmpiricalPairAction &action, const SimulationInfo& simInfo,
  const int norder, const double rmin, const double rmax, const int ngpts,
  const int nImages, int exLevel) 
  : PairAction(s1,s2,action,simInfo,norder,rmin,rmax,ngpts,false,exLevel),
    nImages(nImages) {
  std :: vector<double> cell(NDIM);
  sphereR=0.;
  for (int i=0; i< NDIM; i++){
    cell[i] =  (*simInfo.getSuperCell())[i];
    sphereR = (cell[i]>sphereR)?cell[i]:sphereR;
  }
  sphereR *=nImages;
  boxImageVecs.resize(0);
  findBoxImageVectors(cell);
  std::cout << "Using Trad Ewald with nImages: " << nImages << std::endl;
} 


double EwaldImagePairAction::getActionDifference(const MultiLevelSampler& sampler,
                                         const int level) {
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads();
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads();
  const SuperCell& cell=sampler.getSuperCell();
  const int nStride=(int)pow(2,level);
  const int nSlice=sectionBeads.getNSlice();
  const IArray& index=sampler.getMovingIndex(); 
  const int nMoving=index.size();
  double deltaAction=0;
  for (int iMoving=0; iMoving<nMoving; ++iMoving) {
    const int i=index(iMoving);
    int jbegin,jend;
    if (i>=ifirst1 && i<ifirst1+npart1) {
      jbegin=ifirst2; jend=ifirst2+npart2;
    } else
    if (i>=ifirst2 && i<ifirst2+npart2) {
      jbegin=ifirst1; jend=ifirst1+npart1;
    } else 
    continue; //Particle not in this interaction.
    for (int j=jbegin; j<jend; ++j) {
      bool isMoving=false; int jMoving=0;
      for (int k=0;k<nMoving;++k) {
        if (j==index(k)) {isMoving=true; jMoving=k; break;}
      }
      if (isMoving && i<=j) continue; //Don't double count moving interactions.
      Vec prevDelta =sectionBeads(i,0);
      prevDelta-=sectionBeads(j,0); cell.pbc(prevDelta);
      Vec prevMovingDelta=prevDelta;
      //sum over images
      for (unsigned int img=0; img<boxImageVecs.size(); img++){
	Vec boxImage;
	for (int l=0; l<NDIM; l++) boxImage[l]=boxImageVecs[img][l];

	double prevR=sqrt(dot(prevDelta+boxImage,prevDelta+boxImage));
	double prevMovingR=prevR;
	for (int islice=nStride; islice<nSlice; islice+=nStride) {
	  // Add action for moving beads.
	  Vec delta=movingBeads(iMoving,islice);
	  delta-=(isMoving)?movingBeads(jMoving,islice):sectionBeads(j,islice);
	  cell.pbc(delta);
	  double r=sqrt(dot(delta+boxImage,delta+boxImage));
	  double q=0.5*(r+prevMovingR);
	  if (level<3 && norder>0) {
	    Vec svec=delta-prevMovingDelta; double s2=dot(svec,svec)/(q*q);
	    if (hasZ) {
	      double z=(r-prevMovingR)/q;
	      deltaAction+=uk0(q,s2,z*z);
	    } else {
	      deltaAction+=uk0(q,s2);
	    }
	  } else { 
	    deltaAction+=u00(q);
	  }
	  prevMovingDelta=delta; prevMovingR=r;
	  // Subtract action for old beads.
	  delta=sectionBeads(i,islice);
	  delta-=sectionBeads(j,islice);
	  cell.pbc(delta);
	  r=sqrt(dot(delta+boxImage,delta+boxImage));
	  q=0.5*(r+prevR);
	  if (level<3 && norder>0) {
	    Vec svec=delta-prevDelta; double s2=dot(svec,svec)/(q*q);
	    if (hasZ) {
	      double z=(r-prevR)/q;
	      deltaAction-=uk0(q,s2,z*z);
	    } else {
	      deltaAction-=uk0(q,s2);
	    }
	  } else {
	    deltaAction-=u00(q);
	  }
	  prevDelta=delta; prevR=r;
	}
      }
    }
  }
  return deltaAction*nStride;
}

double EwaldImagePairAction::getActionDifference(const Paths &paths, 
    const VArray &displacement, int nmoving, const IArray &movingIndex, 
    int iFirstSlice, int iLastSlice) {
  const SuperCell& cell=paths.getSuperCell();
  double deltaAction=0;
  for (int iMoving=0; iMoving<nmoving; ++iMoving) {
    const int i=movingIndex(iMoving);
    int jbegin,jend;
    if (i>=ifirst1 && i<ifirst1+npart1) {
      jbegin=ifirst2; jend=ifirst2+npart2;
    } else
    if (i>=ifirst2 && i<ifirst2+npart2) {
      jbegin=ifirst1; jend=ifirst1+npart1;
    } else 
    continue; //Particle not in this interaction.
    for (int j=jbegin; j<jend; ++j) {
      bool isMoving=false; int jMoving=0;
      for (int k=0;k<nmoving;++k) {
        if (j==movingIndex(k)) {isMoving=true; jMoving=k; break;}
      }
      if (isMoving && i<=j) continue; //Don't double count moving interactions.
      Vec prevDelta =paths(i,iFirstSlice-1);
      prevDelta-=paths(j,iFirstSlice-1); cell.pbc(prevDelta);
      Vec prevMovingDelta=prevDelta+displacement(iMoving);
      if (isMoving) prevMovingDelta-=displacement(jMoving);
      cell.pbc(prevMovingDelta);


      for (unsigned int img=0; img<boxImageVecs.size(); img++){//////////
	Vec boxImage;
	for (int l=0; l<NDIM; l++) boxImage[l]=boxImageVecs[img][l];
	
	double prevR=sqrt(dot(prevDelta+boxImage,prevDelta+boxImage));
	double prevMovingR=sqrt(dot(prevMovingDelta+boxImage,prevMovingDelta+boxImage));
	for (int islice=iFirstSlice; islice<=iLastSlice; islice++) {
	  // Add action for moving beads.
	  Vec delta=paths(i,islice);
	  delta+=displacement(iMoving);
	  delta-=paths(j,islice);
	  if (isMoving) delta-=displacement(jMoving);
	  cell.pbc(delta);
	  double r=sqrt(dot(delta+boxImage,delta+boxImage));
	  double q=0.5*(r+prevMovingR);
	  
          Vec svec=delta-prevMovingDelta; double s2=dot(svec,svec)/(q*q);
          deltaAction+=uk0(q,s2);
	  
	  prevMovingDelta=delta; prevMovingR=r;
	  // Subtract action for old beads.
	  delta=paths(i,islice);
	  delta-=paths(j,islice);
	  cell.pbc(delta);
	  r=sqrt(dot(delta+boxImage,delta+boxImage));
	  q=0.5*(r+prevR);
	  
          svec=delta-prevDelta;  s2=dot(svec,svec)/(q*q);
          deltaAction-=uk0(q,s2);
	  
	  prevDelta=delta; prevR=r;
      }
      }
    }
  }
  return deltaAction;
}

double EwaldImagePairAction::getTotalAction(const Paths& paths, int level) const {
  return 0;
}

void EwaldImagePairAction::getBeadAction(const Paths& paths, int ipart, int islice,
    double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const {

  u=utau=ulambda=0; fm=0.; fp=0.;
  int jbegin,jend;
  if (ipart>=ifirst1 && ipart<ifirst1+npart1) {
    jbegin=ifirst2; jend=ifirst2+npart2;
  } else
  if (ipart>=ifirst2 && ipart<ifirst2+npart2) {
    jbegin=ifirst1; jend=ifirst1+npart1;
  } else 
  return; //Particle not in this interaction.


      for (unsigned int img=0; img<boxImageVecs.size(); img++){//////////
	Vec boxImage;
	for (int l=0; l<NDIM; l++) boxImage[l]=boxImageVecs[img][l];


	for (int j=jbegin; j<jend; ++j) {
	  if (ipart==j) continue;
	  Vec delta=paths(ipart,islice);
	  delta-=paths(j,islice);
	  paths.getSuperCell().pbc(delta);
	  double r=sqrt(dot(delta+boxImage,delta+boxImage));
	  Vec prevDelta=paths(ipart,islice,-1);
	  prevDelta-=paths(j,islice,-1);
	  paths.getSuperCell().pbc(prevDelta);
	  double prevR=sqrt(dot(prevDelta+boxImage,prevDelta+boxImage));
	  double q=0.5*(r+prevR);
	  Vec svec=delta-prevDelta; double s2=dot(svec,svec)/(q*q);
	  double v,vtau,vq,vs2,vz2,z;
	  if (hasZ) { 
	    z=(r-prevR)/q;
	    uk0CalcDerivatives(q,s2,z*z,v,vtau,vq,vs2,vz2);
	  } else {
	    uk0CalcDerivatives(q,s2,v,vtau,vq,vs2);
	    vz2=0.; z=0;
	  }
	  u += 0.5*v;
	  utau += 0.5*vtau;
	  
	  
	  fm -= vq*(delta+boxImage)/(2*r) + vs2*(2*svec/(q*q) - s2*(delta+boxImage)/(q*r))
	    +vz2*z*(delta+boxImage)*(2-z)/(q*r);
	  // And force contribution from next slice.
	  Vec nextDelta=paths(ipart,islice,+1);
	  nextDelta-=paths(j,islice,+1);
	  paths.getSuperCell().pbc(nextDelta);
	  double nextR=sqrt(dot(nextDelta+boxImage,nextDelta+boxImage));
	  q=0.5*(r+nextR);
	  svec=delta-nextDelta; s2=dot(svec,svec)/(q*q);
	  if (hasZ) { 
	    z=(r-nextR)/q;
	    uk0CalcDerivatives(q,s2,z*z,v,vtau,vq,vs2,vz2);
	  } else {
	    uk0CalcDerivatives(q,s2,v,vtau,vq,vs2);
	    vz2=0.; z=0.;
	  }
	  fp -= vq*(delta+boxImage)/(2*r) + vs2*(2*svec/(q*q) - s2*(delta+boxImage)/(q*r))
	    +vz2*z*(delta+boxImage)*(2-z)/(q*r);
	}
      }
}



void EwaldImagePairAction::findBoxImageVectors(std :: vector<double> a) { 
 // 3D case first... to be extended to 2D and 1D soon...after testing 3D case
  std :: vector<std :: vector<double> > vertices(8); 
  for (unsigned int i=0; i< vertices.size(); i++){
    vertices[i].resize(NDIM);
  }
  vertices[0][0]=-a[0]/2;  vertices[0][1]=-a[1]/2;  vertices[0][2]=-a[2]/2;
  vertices[1][0]=a[0]/2;  vertices[1][1]=-a[1]/2;  vertices[1][2]=-a[2]/2;
  vertices[2][0]=a[0]/2;  vertices[2][1]=a[1]/2;  vertices[2][2]=-a[2]/2;
  vertices[3][0]=-a[0]/2;  vertices[3][1]=a[1]/2;  vertices[3][2]=-a[2]/2;
  vertices[4][0]=-a[0]/2;  vertices[4][1]=-a[1]/2;  vertices[4][2]=a[2]/2;
  vertices[5][0]=a[0]/2;  vertices[5][1]=-a[1]/2;  vertices[5][2]=a[2]/2;
  vertices[6][0]=a[0]/2;  vertices[6][1]=a[1]/2;  vertices[6][2]=a[2]/2;
  vertices[7][0]=-a[0]/2;  vertices[7][1]=a[1]/2;  vertices[7][2]=a[2]/2;

  // Check that all the box images are inside the sphere of radius nImages*max(cell side)
  // and save the location of these boxes (center of box).
  std :: vector<double> L(NDIM);
  for (int nx=-nImages; nx<=nImages; nx++){
    for (int ny=-nImages; ny<=nImages; ny++){
      for (int nz=-nImages; nz<=nImages; nz++){
	
	L[0]=nx*a[0]; L[1]=ny*a[1]; L[2]=nz*a[2];
	int flag=1;
	for (int v=0;v<8;v++){
	  double tmpR=0;
	  for (int i=0; i< NDIM;i++) {
	    tmpR+=(L[i]+vertices[v][i])*(L[i]+vertices[v][i]);
	  }
	  if (sqrt(tmpR) > sphereR) {
	    flag=0;
	    break;
	  }
	}
	if (flag==1) {
	  boxImageVecs.push_back(L);
	}
	
      }
    }
  }  
}
