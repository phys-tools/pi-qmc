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
#include "AnisotropicNodes.h"
#include "sampler/DoubleMLSampler.h"
#include "SimulationInfo.h"
#include "sampler/DoubleSectionChooser.h"
#include "Species.h"
#include "Paths.h"
#include "Beads.h"
#include "SuperCell.h"
#include <iostream>
#include <cstdlib>
#include <blitz/tinyvec-et.h>

#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
extern "C" void DGETRF_F77(const int*, const int*, double*, const int*,
                           const int*, int*);
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
extern "C" void DGETRI_F77(const int*, double*, const int*, const int*,
                           double*, const int*, int*);

AnisotropicNodes::AnisotropicNodes(const SimulationInfo& simInfo,
  const Species& species, const double temperature,
  const int maxlevel, const int maxmovers)
  : temperature(temperature), tau(simInfo.getTau()),
    mass(*species.anMass),sqrtmass(sqrt(mass)),
    npart(species.count),ifirst(species.ifirst),
    r1(npart), r2(npart),
    matrix((int)(pow(2,maxlevel+1)+0.1)+1),
    newMatrix((int)(pow(2,maxlevel+1)+0.1)+1),
    newColumn((int)(pow(2,maxlevel+1)+0.1)+1),
    dist((int)(pow(2,maxlevel+1)+0.1)+1,npart),
    newDist((int)(pow(2,maxlevel+1)+0.1)+1,npart),
    det((int)(pow(2,maxlevel+1)+0.1)+1),
    newDet((int)(pow(2,maxlevel+1)+0.1)+1),
    ipiv(npart),lwork(npart*npart),
    work(lwork),
    cell(*simInfo.getSuperCell()), 
    notMySpecies(true) {
  for (int idim=0; idim<NDIM; ++idim) {
    pg[idim]=new PeriodicGaussian(mass[idim]*temperature,cell.a[idim],1000);
  }
  for (unsigned int i=0; i<matrix.size(); ++i) { 
    matrix[i]=new Matrix(npart,npart,ColMajor());
    newMatrix[i]=new Matrix(npart,npart,ColMajor());
    newColumn[i]=new Matrix(npart,maxmovers,ColMajor());
  }
  std::cout << "AnisotropicNodes with temperature = "
            << temperature << std::endl;
}

AnisotropicNodes::~AnisotropicNodes() {
  for (int i=0; i<NDIM; ++i) delete pg[i];
  for (unsigned int i=0; i<matrix.size(); ++i) {
    delete matrix[i];
    delete newMatrix[i];
    delete newColumn[i];
  }
}

void AnisotropicNodes::initialize(const DoubleSectionChooser& chooser) {
  sectionchooser=&chooser;
  // Get ready to move paths.
  const Beads<NDIM>& sectionBeads=chooser.getBeads(1);
  const Beads<NDIM>& otherBeads=chooser.getBeads(2);
  const int nSlice=sectionBeads.getNSlice();
  for (int islice=0; islice<nSlice; ++islice) {
    for (int i=0; i<npart; ++i) r2(i)=otherBeads(i+ifirst,islice);
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads(i+ifirst,islice);
    det(islice)=evaluate(r1,r2,*matrix[islice]);
    if (det(islice)<0.0) std::cout << "ERROR - det<0" << std::endl;
    evaluateDistances(islice);
  } 
} 


double AnisotropicNodes::getActionDifference(
    const DoubleMLSampler& sampler, int level) {
  // Get ready to move paths.
  double deltaAction=0;
  movingIndex=&sampler.getMovingIndex(1); 
  const IArray index=(*movingIndex); 
  if (index(0)<ifirst || index(0)>=ifirst+npart) {
    notMySpecies=true; 
    return 0;
  }
  notMySpecies=false;
  const Beads<NDIM>& sectionBeads=sampler.getSectionBeads(1);
  const Beads<NDIM>& otherBeads=sampler.getSectionBeads(2);
  const Beads<NDIM>& movingBeads=sampler.getMovingBeads(1);
  const int nSlice=sectionBeads.getNSlice();
  const int nMoving=index.size();
  const int nStride=(int)pow(2,level+1);
  // First check for node crossing.
  for (int islice=nStride/2; islice<nSlice; islice+=nStride) {
    for (int i=0; i<npart; ++i) r2(i)=otherBeads(i+ifirst,islice);
    for (int i=0; i<npart; ++i) r1(i)=sectionBeads(i+ifirst,islice);
    for (int i=0; i<nMoving; ++i) r1(index(i)-ifirst)=movingBeads(i,islice);
    newDet(islice)=evaluateUpdate(r1,r2,islice);
    if (newDet(islice)<0) return deltaAction=2e100;
  } 
  // Calculate the nodal action if level=0;
  if (level==0) {
    const double twomovert=2.0/tau; //move factor of mass in the distance^2
    double d0(0),newd0(0),d1(0),newd1(0);
    for (int islice=0;islice<nSlice; ++islice) { 
      *newMatrix[islice]=*matrix[islice];
      d1=0; for (int i=0;i<npart;++i) d1+=1./(dist(islice,i)*dist(islice,i));
      newd1=d1=sqrt(1./d1);
      if (islice>0) {
        for (int i=0; i<npart; ++i) r2(i)=otherBeads(i+ifirst,islice);
        for (int i=0; i<npart; ++i) r1(i)=sectionBeads(i+ifirst,islice);
        for (int i=0; i<nMoving; ++i) r1(index(i)-ifirst)=movingBeads(i,islice);
        newd1=evaluateGrad(islice);
        deltaAction+=-log(1-exp(-twomovert*newd1*newd0))
                     +log(1-exp(-twomovert*d1*d0));
	//std::cout << newd1 << " " << newd0 << " " << d1 << " " << d0 << std::endl;
      }
      d0=d1; newd0=newd1; 
    }
  }
  return deltaAction;
}

double AnisotropicNodes::evaluateGrad(const int islice) {
  // First update the slater matrix.
  const IArray index=(*movingIndex); 
  const int nmoving=index.size();
  Matrix& mat(*newMatrix[islice]);
  Matrix& col(*newColumn[islice]);
  //VMatrix& colGrad(*newColGrad[islice]);
  for (int kmoving=0; kmoving<nmoving; ++kmoving) {
    int kpart=index(kmoving)-ifirst;
    if (kpart<0 || kpart>=npart) break;
    double detratio=0;
    for (int jpart=0; jpart<npart; ++jpart) {
      detratio+=mat(kpart,jpart)*col(jpart,kmoving);
    }
    for (int jpart=0; jpart<npart; ++jpart) mat(kpart,jpart)/=detratio;  
    for (int ipart=0; ipart<npart; ++ipart) {
      if (ipart!=kpart) {
        double sum=0;
        for (int lpart=0; lpart<npart; ++lpart) {
          sum+=mat(ipart,lpart)*col(lpart,kmoving); 
        }
        for (int jpart=0; jpart<npart; ++jpart) {
          mat(ipart,jpart)-=mat(kpart,jpart)*sum; 
        }
      }
    }
  }
  /// Next evaluate the gradients.
  double d2inv=0;
  for(int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0;
    for(int ipart=0; ipart<npart; ++ipart) {
      Vec delta=r1(jpart)-r2(ipart);
      cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=pg[i]->grad(fabs(delta[i]))*sqrtmass[i];
        if (delta[i]<0) grad[i]=-grad[i];
      }
      logGrad+=mat(jpart,ipart)*grad;
    }
    newDist(islice,jpart)=1.0/sqrt(dot(logGrad,logGrad));
    d2inv+=1./(newDist(islice,jpart)*newDist(islice,jpart));
  }
  return sqrt(1./d2inv);
} 

void AnisotropicNodes::evaluateDistances(const int islice) {
  const Matrix& mat(*matrix[islice]);
  for(int jpart=0; jpart<npart; ++jpart) {
    Vec logGrad=0.0;
    for(int ipart=0; ipart<npart; ++ipart) {
      Vec delta=r1(jpart)-r2(ipart);
      cell.pbc(delta);
      Vec grad;
      for (int i=0; i<NDIM; ++i) {
        grad[i]=pg[i]->grad(fabs(delta[i]))*sqrtmass[i];
        if (delta[i]<0) grad[i]=-grad[i];
      }
      logGrad+=mat(jpart,ipart)*grad;
    }
    dist(islice,jpart)=1.0/sqrt(dot(logGrad,logGrad));
  }
} 

double AnisotropicNodes::evaluate(const VArray& r1, const VArray& r2, 
                                   Matrix& mat) {
  mat=0;
  for(int jpart=0; jpart<npart; ++jpart) {
    for(int ipart=0; ipart<npart; ++ipart) {
      Vec delta=r1(jpart)-r2(ipart);
      cell.pbc(delta);
      double ear2=1;
      for (int i=0; i<NDIM; ++i) ear2*=(*pg[i])(fabs(delta[i]));
      mat(ipart,jpart)=ear2;
    }
  }
  // Calculate determinant and inverse.
  int info=0; 
  //LU decomposition
  DGETRF_F77(&npart,&npart,mat.data(),&npart,ipiv.data(),&info);
  if (info!=0) std::cout << "BAD RETURN FROM ZGETRF!!!!" << std::endl;
  double det = 1;
  for (int i=0; i<npart; ++i) det*= mat(i,i) * ((i+1==ipiv(i))?1:-1);
  //Find inverse.
  DGETRI_F77(&npart,mat.data(),&npart,ipiv.data(),work.data(),&lwork,&info);
  return det;
} 

double AnisotropicNodes::evaluateUpdate(const VArray& r1, const VArray& r2, 
    const int islice) {
  const IArray index=(*movingIndex); 
  int nmoving=index.size();
  Matrix& mat(*matrix[islice]);
  Matrix& col(*newColumn[islice]);
  Matrix detMat(nmoving,nmoving,ColMajor()); detMat=0;
  for (int kmoving=0; kmoving<nmoving; ++kmoving) {
    int kpart=index(kmoving)-ifirst;
    if (kpart<0 || kpart>=npart) {detMat(kmoving,kmoving)=1.0; break;}
    for (int jpart=0; jpart<npart; ++jpart) {
      Vec delta=r1(kpart)-r2(jpart);
      cell.pbc(delta);
      double ear2=1; for (int i=0;i<NDIM;++i) ear2*=(*pg[i])(fabs(delta[i]));
      col(jpart,kmoving)=ear2;
      for (int imoving=0; imoving<nmoving; ++imoving) {
        int ipart=index(imoving)-ifirst;
        if (ipart<0 || ipart>=npart) break;
        detMat(imoving,kmoving) += mat(ipart,jpart)*col(jpart,kmoving);
      }
    }
  }
  // Calculate determinant and inverse.
  int info=0;
  dgetrf_(&nmoving,&nmoving,detMat.data(),&nmoving,ipiv.data(),&info);
  if (info!=0) std::cout << "BAD RETURN from ZGETRF!!!" << std::endl;
  double newdet=det(islice);
  for (int i=0; i<nmoving; ++i) newdet*= detMat(i,i)*((i+1==ipiv(i))?1:-1);
  return newdet;
}

double AnisotropicNodes::getTotalAction(const Paths&, const int level) const {
  return 0;
}

void AnisotropicNodes::getBeadAction(const Paths&, int ipart, int islice,
           double& u, double& utau, double& ulambda, Vec& fp, Vec& fm) const {
}

void AnisotropicNodes::acceptLastMove() {
  const IArray index=(*movingIndex); 
  if (notMySpecies) return;
  for (int i=1; i<det.size(); ++ i) {
    Matrix* mat=matrix[i]; matrix[i]=newMatrix[i]; newMatrix[i]=mat;
    double d=det(i); det(i)=newDet(i); newDet(i)=d;
    for (int j=0; j<npart; ++j) dist(i,j)=newDist(i,j);
  }
}
